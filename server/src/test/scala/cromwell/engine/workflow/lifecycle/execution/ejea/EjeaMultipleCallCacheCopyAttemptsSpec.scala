package cromwell.engine.workflow.lifecycle.execution.ejea

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

class EjeaMultipleCallCacheCopyAttemptsSpec
  extends EngineJobExecutionActorSpec
    with HasJobSuccessResponse
    with HasCopyFailureResponses
    with HasJobFailureResponses
    with CanExpectJobStoreWrites
    with CanExpectFetchCachedResults {

  override implicit val stateUnderTest = BackendIsCopyingCachedOutputs
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
      val response = copyAttemptFailedResponse(copyAttemptNumber)

      ejea ! response
      helper.ejhaProbe.expectNoMessage()

      eventually {
        ejea.stateName should be(RunningJob)
      }
    }

    def copyAttemptSucceedsAndEjeaReacts() = {
      ejea ! successResponse
      expectJobStoreWrite(expectedData = SucceededResponseData(successResponse, None))
    }

    "keep waiting for success up to the configured max-failed-copy-attempts limit" in {
      val maxFailedCopyAttempts = 100
      ejea = buildEjea(maxFailedCopyAttempts)

      // First: A long series of copy failure:
      val initialCopyFailures = maxFailedCopyAttempts - 1
      0.until(initialCopyFailures). foreach { currentCopyAttemptNumber =>
        ejhaSendsHitIdToEjeaAndEjeaReacts(currentCopyAttemptNumber)
        fetchCachedResultsActorSendsResultSetToEjeaAndEjeaReacts(currentCopyAttemptNumber)
        copyAttemptFailsAndEjeaLooksForNextHit(becauseBlacklisted = false, currentCopyAttemptNumber)
      }

      // Then: A success:
      val currentCopyAttemptNumber = initialCopyFailures // because the initial 0.until(...) is non-inclusive of the argument
      ejhaSendsHitIdToEjeaAndEjeaReacts(copyAttemptNumber = currentCopyAttemptNumber)
      fetchCachedResultsActorSendsResultSetToEjeaAndEjeaReacts(currentCopyAttemptNumber)
      copyAttemptSucceedsAndEjeaReacts()
    }

    "fail fast after the configured max-failed-copy-attempts limit is hit" in {
      val maxFailedCopyAttempts = 100
      ejea = buildEjea(maxFailedCopyAttempts)

      // First: A long series of copy failure:
      val initialCopyFailures = maxFailedCopyAttempts - 1
      0.until(initialCopyFailures). foreach { currentCopyAttemptNumber =>
        ejhaSendsHitIdToEjeaAndEjeaReacts(currentCopyAttemptNumber)
        fetchCachedResultsActorSendsResultSetToEjeaAndEjeaReacts(currentCopyAttemptNumber)
        copyAttemptFailsAndEjeaLooksForNextHit(becauseBlacklisted = false, currentCopyAttemptNumber)
      }

      // Then: Another failure:
      val currentCopyAttemptNumber = initialCopyFailures // because the initial 0.until(...) is non-inclusive of the argument
      ejhaSendsHitIdToEjeaAndEjeaReacts(copyAttemptNumber = currentCopyAttemptNumber)
      fetchCachedResultsActorSendsResultSetToEjeaAndEjeaReacts(currentCopyAttemptNumber)
      copyAttemptFailsAndEjeaGivesUp(currentCopyAttemptNumber)
    }

    "disregard any number of blacklist caused-failures on the road to the max-failed-copy-attempts limit" in {
      val maxFailedCopyAttempts = 100
      ejea = buildEjea(maxFailedCopyAttempts)

      // First: A long series of (genuine) copy failures:
      val initialCopyFailures = maxFailedCopyAttempts - 1
      0.until(initialCopyFailures). foreach { currentCopyAttemptNumber =>
        ejhaSendsHitIdToEjeaAndEjeaReacts(currentCopyAttemptNumber)
        fetchCachedResultsActorSendsResultSetToEjeaAndEjeaReacts(currentCopyAttemptNumber)
        copyAttemptFailsAndEjeaLooksForNextHit(becauseBlacklisted = false, currentCopyAttemptNumber)
      }

      // Second: An even longer series of (blacklist) copy failures:
      val blacklistCopyFailures = maxFailedCopyAttempts + 2
      initialCopyFailures.until(initialCopyFailures + blacklistCopyFailures). foreach { currentCopyAttemptNumber =>
        ejhaSendsHitIdToEjeaAndEjeaReacts(currentCopyAttemptNumber)
        fetchCachedResultsActorSendsResultSetToEjeaAndEjeaReacts(currentCopyAttemptNumber)
        copyAttemptFailsAndEjeaLooksForNextHit(becauseBlacklisted = true, currentCopyAttemptNumber)
      }

      // Then: Another (genuine) failure:
      val currentCopyAttemptNumber = initialCopyFailures + blacklistCopyFailures
      ejhaSendsHitIdToEjeaAndEjeaReacts(copyAttemptNumber = currentCopyAttemptNumber)
      fetchCachedResultsActorSendsResultSetToEjeaAndEjeaReacts(currentCopyAttemptNumber)
      copyAttemptFailsAndEjeaGivesUp(currentCopyAttemptNumber)
    }


  }

  def buildEjea(maxFailedCopyAttempts: Int) = helper
    .buildEJEA(
      restarting = false,
      callCachingMode = CallCachingActivity(ReadCache, CallCachingOptions(invalidateBadCacheResults = false)),
      callCachingMaxFailedCopyAttempts = maxFailedCopyAttempts)
    .setStateInline(state = CheckingCallCache, data = ResponsePendingData(
      jobDescriptor = helper.backendJobDescriptor,
      bjeaProps = helper.bjeaProps,
      ejha = Some(helper.ejhaProbe.ref)))

}
