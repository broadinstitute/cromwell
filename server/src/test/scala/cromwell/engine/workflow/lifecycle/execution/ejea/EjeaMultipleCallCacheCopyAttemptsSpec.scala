package cromwell.engine.workflow.lifecycle.execution.ejea

import com.typesafe.scalalogging.StrictLogging
import cromwell.backend.BackendCacheHitCopyingActor.CopyOutputsCommand
import cromwell.core.WorkflowId
import cromwell.core.callcaching._
import cromwell.core.simpleton.WomValueSimpleton
import cromwell.engine.workflow.lifecycle.execution.callcaching.CallCacheReadingJobActor.NextHit
import cromwell.engine.workflow.lifecycle.execution.callcaching.EngineJobHashingActor.CacheHit
import cromwell.engine.workflow.lifecycle.execution.callcaching.FetchCachedResultsActor.CachedOutputLookupSucceeded
import cromwell.engine.workflow.lifecycle.execution.ejea.EngineJobExecutionActorSpec._
import cromwell.engine.workflow.lifecycle.execution.job.EngineJobExecutionActor._
import cromwell.services.CallCaching.CallCachingEntryId
import wom.values.WomString
import EjeaMultipleCallCacheCopyAttemptsSpec.bt140Debug
import akka.testkit.TestFSMRef

class EjeaMultipleCallCacheCopyAttemptsSpec
  extends EngineJobExecutionActorSpec
    with HasJobSuccessResponse
    with HasCopyFailureResponses
    with HasJobFailureResponses
    with CanExpectJobStoreWrites
    with CanExpectFetchCachedResults {

  override implicit val stateUnderTest: EngineJobExecutionActorState = BackendIsCopyingCachedOutputs
  override val allowMultipleCacheCycles: Boolean = true

  "An EJEA attempting to call cache copy" should {

    // Arbitrary.
    // When we attempt the nth copy attempt, we'll say that the cache entry ID is 'n' plus this offset.
    // Just makes sure that we're treating the copy attempt and the hit ID as separate numbers.
    def cacheEntryIdFromCopyAttempt(attempt: Int) = CallCachingEntryId(75 + attempt)

    def ejhaSendsHitIdToEjeaAndEjeaReacts(copyAttemptNumber: Int) = {
      val callCachingEntryId = cacheEntryIdFromCopyAttempt(copyAttemptNumber)
      helper.ejhaProbe.send(ejea, CacheHit(callCachingEntryId))
      expectFetchCachedResultsActor(callCachingEntryId)
      eventually {
        ejea.stateName should be(FetchingCachedOutputsFromDatabase)
        ejea.stateData.asInstanceOf[ResponsePendingData].ejeaCacheHit.get.hitNumber should be(copyAttemptNumber)
      }
    }

    def fetchCachedResultsActorSendsResultSetToEjeaAndEjeaReacts(copyAttemptNumber: Int) = {
      val callCachingEntryId = cacheEntryIdFromCopyAttempt(copyAttemptNumber)
      val cachedSimpletons = Seq(WomValueSimpleton("a", WomString("hullo")), WomValueSimpleton("b", WomString("cheerio")))
      val detritusMap = Map("stdout" -> "//somePath")
      val cachedReturnCode = Some(17)
      val sourceCacheDetails = s"${WorkflowId.randomId()}:call-someTask:1"
      ejea ! CachedOutputLookupSucceeded(cachedSimpletons, detritusMap, cachedReturnCode, callCachingEntryId, sourceCacheDetails)
      helper.callCacheHitCopyingProbe.expectMsg(CopyOutputsCommand(cachedSimpletons, detritusMap, callCachingEntryId, cachedReturnCode))
      eventually {
        ejea.stateName should be(BackendIsCopyingCachedOutputs)
      }
    }

    def copyAttemptFailsAndEjeaLooksForNextHit(becauseBlacklisted: Boolean, copyAttemptNumber: Int) = {
      val response = if (becauseBlacklisted) cacheHitBlacklistedResponse(copyAttemptNumber) else copyAttemptFailedResponse(copyAttemptNumber)

      ejea ! response
      helper.ejhaProbe.expectMsg(NextHit)
      eventually {
        ejea.stateName should be(CheckingCallCache)
      }
    }

    def copyAttemptFailsAndEjeaGivesUp(copyAttemptNumber: Int) = {
      bt140Debug("'copyAttemptFailsAndEjeaGivesUp' starting")
      val response = copyAttemptFailedResponse(copyAttemptNumber)

      ejea ! response
      bt140Debug("'copyAttemptFailsAndEjeaGivesUp' expecting no message")
      helper.ejhaProbe.expectNoMessage()

      bt140Debug("'copyAttemptFailsAndEjeaGivesUp' waiting for 'RunningJob'")
      val result = eventually {
        ejea.stateName should be(RunningJob)
      }
      bt140Debug("'copyAttemptFailsAndEjeaGivesUp' done")
      result
    }

    def copyAttemptSucceedsAndEjeaReacts(): Unit = {
      ejea ! successResponse
      expectJobStoreWrite(expectedData = SucceededResponseData(successResponse, None))
    }

    "keep waiting for success up to the configured max-failed-copy-attempts limit" in {
      bt140Debug("'keep waiting' starting")
      val maxFailedCopyAttempts = 100
      ejea = buildEjea(maxFailedCopyAttempts)

      // First: A long series of copy failure:
      bt140Debug("'keep waiting' running 'series of copy failures'")
      val initialCopyFailures = maxFailedCopyAttempts - 1
      0.until(initialCopyFailures). foreach { currentCopyAttemptNumber =>
        ejhaSendsHitIdToEjeaAndEjeaReacts(currentCopyAttemptNumber)
        fetchCachedResultsActorSendsResultSetToEjeaAndEjeaReacts(currentCopyAttemptNumber)
        copyAttemptFailsAndEjeaLooksForNextHit(becauseBlacklisted = false, currentCopyAttemptNumber)
      }

      // Then: A success:
      val currentCopyAttemptNumber = initialCopyFailures // because the initial 0.until(...) is non-inclusive of the argument
      bt140Debug("'keep waiting' running 'ejhaSendsHitIdToEjeaAndEjeaReacts'")
      ejhaSendsHitIdToEjeaAndEjeaReacts(copyAttemptNumber = currentCopyAttemptNumber)
      bt140Debug("'keep waiting' running 'fetchCachedResultsActorSendsResultSetToEjeaAndEjeaReacts'")
      fetchCachedResultsActorSendsResultSetToEjeaAndEjeaReacts(currentCopyAttemptNumber)
      bt140Debug("'keep waiting' running 'copyAttemptSucceedsAndEjeaReacts'")
      copyAttemptSucceedsAndEjeaReacts()
      bt140Debug("'keep waiting' done")
    }

    "fail fast after the configured max-failed-copy-attempts limit is hit" in {
      bt140Debug("'fail fast' starting")
      val maxFailedCopyAttempts = 100
      ejea = buildEjea(maxFailedCopyAttempts)

      // First: A long series of copy failure:
      bt140Debug("'fail fast' running 'longer series of (genuine) copy failures'")
      val initialCopyFailures = maxFailedCopyAttempts - 1
      0.until(initialCopyFailures). foreach { currentCopyAttemptNumber =>
        ejhaSendsHitIdToEjeaAndEjeaReacts(currentCopyAttemptNumber)
        fetchCachedResultsActorSendsResultSetToEjeaAndEjeaReacts(currentCopyAttemptNumber)
        copyAttemptFailsAndEjeaLooksForNextHit(becauseBlacklisted = false, currentCopyAttemptNumber)
      }

      // Then: Another failure:
      bt140Debug("'fail fast' running 'ejhaSendsHitIdToEjeaAndEjeaReacts'")
      val currentCopyAttemptNumber = initialCopyFailures // because the initial 0.until(...) is non-inclusive of the argument
      ejhaSendsHitIdToEjeaAndEjeaReacts(copyAttemptNumber = currentCopyAttemptNumber)
      bt140Debug("'fail fast' running 'fetchCachedResultsActorSendsResultSetToEjeaAndEjeaReacts'")
      fetchCachedResultsActorSendsResultSetToEjeaAndEjeaReacts(currentCopyAttemptNumber)
      bt140Debug("'fail fast' running 'copyAttemptFailsAndEjeaGivesUp'")
      copyAttemptFailsAndEjeaGivesUp(currentCopyAttemptNumber)
      bt140Debug("'fail fast' done")
    }

    "disregard any number of blacklist caused-failures on the road to the max-failed-copy-attempts limit" in {
      bt140Debug("'disregard' starting")
      val maxFailedCopyAttempts = 100
      ejea = buildEjea(maxFailedCopyAttempts)

      // First: A long series of (genuine) copy failures:
      bt140Debug("'disregard' running 'longer series of (genuine) copy failures'")
      val initialCopyFailures = maxFailedCopyAttempts - 1
      0.until(initialCopyFailures). foreach { currentCopyAttemptNumber =>
        ejhaSendsHitIdToEjeaAndEjeaReacts(currentCopyAttemptNumber)
        fetchCachedResultsActorSendsResultSetToEjeaAndEjeaReacts(currentCopyAttemptNumber)
        copyAttemptFailsAndEjeaLooksForNextHit(becauseBlacklisted = false, currentCopyAttemptNumber)
      }

      // Second: An even longer series of (blacklist) copy failures:
      bt140Debug("'disregard' running 'longer series of (exclude list) copy failures'")
      val blacklistCopyFailures = maxFailedCopyAttempts + 2
      initialCopyFailures.until(initialCopyFailures + blacklistCopyFailures). foreach { currentCopyAttemptNumber =>
        ejhaSendsHitIdToEjeaAndEjeaReacts(currentCopyAttemptNumber)
        fetchCachedResultsActorSendsResultSetToEjeaAndEjeaReacts(currentCopyAttemptNumber)
        copyAttemptFailsAndEjeaLooksForNextHit(becauseBlacklisted = true, currentCopyAttemptNumber)
      }

      // Then: Another (genuine) failure:
      val currentCopyAttemptNumber = initialCopyFailures + blacklistCopyFailures
      bt140Debug("'disregard' running 'ejhaSendsHitIdToEjeaAndEjeaReacts'")
      ejhaSendsHitIdToEjeaAndEjeaReacts(copyAttemptNumber = currentCopyAttemptNumber)
      bt140Debug("'disregard' running 'fetchCachedResultsActorSendsResultSetToEjeaAndEjeaReacts'")
      fetchCachedResultsActorSendsResultSetToEjeaAndEjeaReacts(currentCopyAttemptNumber)
      bt140Debug("'disregard' running 'copyAttemptFailsAndEjeaGivesUp'")
      copyAttemptFailsAndEjeaGivesUp(currentCopyAttemptNumber)
      bt140Debug("'disregard' done")
    }


  }

  def buildEjea(maxFailedCopyAttempts: Int): TestFSMRef[EngineJobExecutionActorState, EJEAData, MockEjea] = helper
    .buildEJEA(
      restarting = false,
      callCachingMode = CallCachingActivity(ReadCache, CallCachingOptions(invalidateBadCacheResults = false)),
      callCachingMaxFailedCopyAttempts = maxFailedCopyAttempts)
    .setStateInline(state = CheckingCallCache, data = ResponsePendingData(
      jobDescriptor = helper.backendJobDescriptor,
      bjeaProps = helper.bjeaProps,
      ejha = Some(helper.ejhaProbe.ref)))

}

object EjeaMultipleCallCacheCopyAttemptsSpec extends StrictLogging {
  /*
  During the BT-140 investigation the `EjeaMultipleCallCacheCopyAttemptsSpec` was seen not responding in one of the test
  runs. This could have been a GC pause, Akka TestKit not timing out as expected, or any number of other reasons. For
  now, we're adding a little bit more debug information to this spec just in case this happens again in the future. If
  at any point someone decides that this is not actually helping, feel free to delete all calls to this method and the
  debug method itself.
   */
  private def bt140Debug(message: String): Unit = {
    logger.info("BT-140 debug: " + message)
  }
}
