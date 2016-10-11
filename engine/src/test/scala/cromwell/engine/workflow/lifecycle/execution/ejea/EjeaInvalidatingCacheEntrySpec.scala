package cromwell.engine.workflow.lifecycle.execution.ejea

import cats.data.NonEmptyList
import cromwell.core.callcaching.{CallCachingActivity, ReadCache}
import cromwell.engine.workflow.lifecycle.execution.EngineJobExecutionActor._
import cromwell.engine.workflow.lifecycle.execution.callcaching.EngineJobHashingActor.CacheHit
import cromwell.engine.workflow.lifecycle.execution.callcaching.{CallCacheInvalidatedFailure, CallCacheInvalidatedSuccess, MetaInfoId}
import cromwell.engine.workflow.lifecycle.execution.ejea.EngineJobExecutionActorSpec._

class EjeaInvalidatingCacheEntrySpec extends EngineJobExecutionActorSpec with CanExpectFetchCachedResults {

  override implicit val stateUnderTest = InvalidatingCacheEntry

  "An EJEA in InvalidatingCacheEntry state" should {

    val invalidationErrorCause = new Exception("blah")
    val invalidateSuccess = CallCacheInvalidatedSuccess
    val invalidateFailure = CallCacheInvalidatedFailure(invalidationErrorCause)

    val metaInfo31: MetaInfoId = MetaInfoId(31)
    val metaInfo32: MetaInfoId = MetaInfoId(32)
    val cacheHitWithTwoIds = CacheHit(NonEmptyList.of(metaInfo32, metaInfo31))
    val cacheHitWithSingleId = CacheHit(NonEmptyList.of(metaInfo31))
    val noCacheHit = None

    case class CacheHitAndNextExpected(name: String,
                                       cacheHitData: CacheHit,
                                       nextCacheHitToTry: MetaInfoId,
                                       nextExpectedCacheHitData: Option[CacheHit])

    val cacheHitCombinations = List(
      CacheHitAndNextExpected("cache hit with single Id hit", cacheHitWithSingleId, metaInfo31, None),
      CacheHitAndNextExpected("cache hit with two Id hits", cacheHitWithTwoIds, metaInfo32, Option(cacheHitWithSingleId))
    )

    RestartOrExecuteCommandTuples foreach { case RestartOrExecuteCommandTuple(operationName, restarting, expectedMessage) =>
      List(invalidateSuccess, invalidateFailure) foreach { invalidateActorResponse =>
        s"$operationName a job if cache invalidation succeeds and there are no other cache hits to try when invalidate response is $invalidateActorResponse" in {
          ejea = ejeaInvalidatingCacheEntryState(noCacheHit, restarting = restarting)
          // Send the response from the invalidate actor
          ejea ! invalidateActorResponse

          helper.bjeaProbe.expectMsg(awaitTimeout, expectedMessage)
          ejea.stateName should be(RunningJob)
          ejea.stateData should be(ResponsePendingData(helper.backendJobDescriptor, helper. bjeaProps, None, None))
        }
      }
    }

    cacheHitCombinations foreach { case combo @ CacheHitAndNextExpected(hitComboName, initialHitData, nextHitToTry, finalHitData) =>
        List(invalidateSuccess, invalidateFailure) foreach { invalidateActorResponse =>
          s"try the next available hit when there is a $hitComboName when response is $invalidateActorResponse" in {
            ejea = ejeaInvalidatingCacheEntryState(Option(initialHitData), restarting = false)
            // Send the response from the invalidate actor
            ejea ! invalidateActorResponse

            helper.bjeaProbe.expectNoMsg(awaitAlmostNothing)
            expectFetchCachedResultsActor(nextHitToTry)
            ejea.stateName should be(FetchingCachedOutputsFromDatabase)
            ejea.stateData should be(ResponsePendingData(helper.backendJobDescriptor, helper.bjeaProps, None, finalHitData))
          }
        }
    }
  }

  def standardResponsePendingData(hit: Option[CacheHit]) = ResponsePendingData(helper.backendJobDescriptor, helper.bjeaProps, None, hit)
  def ejeaInvalidatingCacheEntryState(hit: Option[CacheHit], restarting: Boolean = false) = helper.buildEJEA(restarting = restarting, callCachingMode = CallCachingActivity(ReadCache)).setStateInline(data = standardResponsePendingData(hit))
}
