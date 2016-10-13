package cromwell.engine.workflow.lifecycle.execution.ejea

import cats.data.NonEmptyList
import cromwell.core.callcaching.{CallCachingActivity, ReadCache}
import cromwell.engine.workflow.lifecycle.execution.EngineJobExecutionActor._
import cromwell.engine.workflow.lifecycle.execution.callcaching.EngineJobHashingActor.CacheHit
import cromwell.engine.workflow.lifecycle.execution.callcaching.{CallCacheInvalidatedFailure, CallCacheInvalidatedSuccess, CallCachingEntryId}
import cromwell.engine.workflow.lifecycle.execution.ejea.EngineJobExecutionActorSpec._

class EjeaInvalidatingCacheEntrySpec extends EngineJobExecutionActorSpec with CanExpectFetchCachedResults {

  override implicit val stateUnderTest = InvalidatingCacheEntry

  "An EJEA in InvalidatingCacheEntry state" should {

    val invalidationErrorCause = new Exception("blah")
    val invalidateSuccess = CallCacheInvalidatedSuccess
    val invalidateFailure = CallCacheInvalidatedFailure(invalidationErrorCause)

    val metaInfo31: CallCachingEntryId = CallCachingEntryId(31)
    val metaInfo32: CallCachingEntryId = CallCachingEntryId(32)
    val cacheHitWithTwoIds = CacheHit(NonEmptyList.of(metaInfo32, metaInfo31))
    val cacheHitWithSingleId = CacheHit(NonEmptyList.of(metaInfo31))

    List(invalidateSuccess, invalidateFailure) foreach { invalidateActorResponse =>
      s"try the next available hit when response is $invalidateActorResponse" in {
        ejea = ejeaInvalidatingCacheEntryState(Option(cacheHitWithTwoIds), restarting = false)
        // Send the response from the invalidate actor
        ejea ! invalidateActorResponse

        helper.bjeaProbe.expectNoMsg(awaitAlmostNothing)
        expectFetchCachedResultsActor(cacheHitWithSingleId.cacheResultIds.head)
        eventually { ejea.stateName should be(FetchingCachedOutputsFromDatabase) }
        ejea.stateData should be(ResponsePendingData(helper.backendJobDescriptor, helper.bjeaProps, None, Option(cacheHitWithSingleId)))
      }

      RestartOrExecuteCommandTuples foreach { case RestartOrExecuteCommandTuple(operationName, restarting, expectedMessage) =>
        s"$operationName a job if cache invalidation succeeds and there are no other cache hits to try when invalidate response is $invalidateActorResponse" in {
          ejea = ejeaInvalidatingCacheEntryState(Option(cacheHitWithSingleId), restarting = restarting)
          // Send the response from the invalidate actor
          ejea ! invalidateActorResponse

          helper.bjeaProbe.expectMsg(awaitTimeout, expectedMessage)
          eventually { ejea.stateName should be(RunningJob) }
          ejea.stateData should be(ResponsePendingData(helper.backendJobDescriptor, helper. bjeaProps, None, None))
        }
      }
    }
  }

  def standardResponsePendingData(hit: Option[CacheHit]) = ResponsePendingData(helper.backendJobDescriptor, helper.bjeaProps, None, hit)
  def ejeaInvalidatingCacheEntryState(hit: Option[CacheHit], restarting: Boolean = false) = helper.buildEJEA(restarting = restarting, callCachingMode = CallCachingActivity(ReadCache)).setStateInline(data = standardResponsePendingData(hit))
}
