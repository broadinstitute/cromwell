package cromwell.engine.workflow.lifecycle.execution.ejea

import cromwell.backend.BackendCacheHitCopyingActor.{CopyingOutputsFailedResponse, CopyAttemptError, BlacklistSkip}
import cromwell.backend.{BackendJobDescriptor, MetricableCacheCopyErrorCategory}
import cromwell.backend.BackendJobExecutionActor._
import cromwell.core.callcaching._
import cromwell.engine.workflow.lifecycle.execution.callcaching.EngineJobHashingActor.{CallCacheHashes, FileHashes}
import cromwell.engine.workflow.lifecycle.execution.job.EngineJobExecutionActor.{EJEAData, SucceededResponseData, UpdatingCallCache, UpdatingJobStore}
import cromwell.jobstore.JobStoreActor.RegisterJobCompleted
import cromwell.jobstore.{JobResultFailure, JobResultSuccess, JobStoreKey}
import cromwell.services.CallCaching.CallCachingEntryId
import cromwell.util.WomMocks
import org.scalatest.concurrent.Eventually
import wom.values.{WomInteger, WomString}

import scala.util.Success
import scala.util.control.NoStackTrace

private[ejea] trait CanValidateJobStoreKey { self: EngineJobExecutionActorSpec =>
  def validateJobStoreKey(jobKey: JobStoreKey) = {
    jobKey.workflowId should be(helper.workflowId)
    jobKey.callFqn should be(helper.jobFqn)
    jobKey.index should be(helper.jobIndex)
    jobKey.attempt should be(helper.jobAttempt)
  }
}

private[ejea] trait CanExpectCacheWrites extends Eventually { self: EngineJobExecutionActorSpec =>
  def expectCacheWrite(expectedResponse: JobSucceededResponse, expectedCallCacheHashes: CallCacheHashes): Unit = {
    eventually { ejea.stateName should be(UpdatingCallCache) }
    ejea.stateData should be(SucceededResponseData(expectedResponse, Some(Success(expectedCallCacheHashes))))
    ()
  }
}

private[ejea] trait CanExpectJobStoreWrites extends CanValidateJobStoreKey { self: EngineJobExecutionActorSpec with HasJobSuccessResponse with HasJobFailureResponses =>

  def expectJobStoreWrite(expectedData: EJEAData): Unit = {
    helper.jobStoreProbe.expectMsgPF(max = awaitTimeout, hint = "Job Store Write") {
      case RegisterJobCompleted(jobKey, JobResultSuccess(returnCode, jobOutputs)) =>
        validateJobStoreKey(jobKey)
        returnCode should be(successRc)
        jobOutputs should be(successOutputs)
        ejea.stateName should be(UpdatingJobStore)
        ejea.stateData should be(expectedData)
    }
    ()
  }

  def expectJobStoreWriteFailed(expectedData: EJEAData, expectedRetryable: Boolean): Unit = {
    helper.jobStoreProbe.expectMsgPF(max = awaitTimeout, hint = "Job Store Write") {
      case RegisterJobCompleted(jobKey, JobResultFailure(returnCode, reason, retryable)) =>
        validateJobStoreKey(jobKey)
        returnCode should be(failedRc)
        reason should be(failureReason)
        retryable shouldBe expectedRetryable
        ejea.stateName should be(UpdatingJobStore)
        ejea.stateData should be(expectedData)
    }
    ()
  }
}

private[ejea] trait CanExpectHashingInitialization extends Eventually { self: EngineJobExecutionActorSpec =>
  def expectHashingActorInitialization(mode: CallCachingMode, jobDescriptor: BackendJobDescriptor): Unit = {
    eventually { helper.jobHashingInitializations.hasExactlyOne should be(true) }
    helper.jobHashingInitializations.checkIt { initialization =>
      initialization._1 should be(jobDescriptor)
      initialization._2 should be(mode)
    }
  }
}

private[ejea] trait CanExpectFetchCachedResults extends Eventually { self: EngineJobExecutionActorSpec =>
  def expectFetchCachedResultsActor(expectedCallCachingEntryId: CallCachingEntryId): Unit = {
    eventually {
      if (allowMultipleCacheCycles) { helper.fetchCachedResultsActorCreations.hasAtLeastOne should be(true) }
      else { helper.fetchCachedResultsActorCreations.hasExactlyOne should be(true) }
    }
    helper.fetchCachedResultsActorCreations.checkLatest {
      case (callCachingEntryId, _) => callCachingEntryId should be(expectedCallCachingEntryId)
      case _ => fail("Incorrect creation of the fetchCachedResultsActor")
    }
  }
}

private[ejea] trait CanExpectCacheInvalidation extends Eventually { self: EngineJobExecutionActorSpec =>
  def expectInvalidateCallCacheActor(expectedCacheId: CallCachingEntryId): Unit = {
    eventually { helper.invalidateCacheActorCreations.hasExactlyOne should be(true) }
    helper.invalidateCacheActorCreations.checkIt { cacheId =>
      cacheId shouldBe expectedCacheId
    }
  }
}

private[ejea] trait HasJobSuccessResponse { self: EngineJobExecutionActorSpec =>
  val successRc = Option(171)
  val successOutputs = WomMocks.mockOutputExpectations(Map("a" -> WomInteger(3), "b" -> WomString("bee")))
  def successResponse = JobSucceededResponse(helper.jobDescriptorKey, successRc, successOutputs, None, Seq.empty, None, resultGenerationMode = RunOnBackend)
}
private[ejea] object HasJobSuccessResponse {
  val SuccessfulCallCacheHashes = CallCacheHashes(
    Set(HashResult(HashKey("whatever you want"), HashValue("whatever you need"))),
    "initialHash",
    Option(FileHashes(
      Set(HashResult(HashKey("whatever file you want"), HashValue("whatever file you need"))),
      "fileHash"
    ))
  )
}

private[ejea] trait HasJobFailureResponses { self: EngineJobExecutionActorSpec =>
  val failedRc = Option(12)
  val failureReason = new Exception("Deliberate failure for test case: job run failed!") with NoStackTrace
  // Need to delay making the response because job descriptors come from the per-test "helper", which is null outside tests!
  def failureRetryableResponse = JobFailedRetryableResponse(helper.jobDescriptorKey, failureReason, failedRc)
  def failureNonRetryableResponse = JobFailedNonRetryableResponse(helper.jobDescriptorKey, failureReason, Option(12))
  def abortedResponse = JobAbortedResponse(helper.jobDescriptorKey)
}

private[ejea] trait HasCopyFailureResponses { self: EngineJobExecutionActorSpec =>
  val copyFailureReason =
    new Exception("Deliberate failure for test case: failed to copy cache outputs!") with NoStackTrace

  // Need to delay making the response because job descriptors come from the per-test "helper", which is null outside tests!
  def copyAttemptFailedResponse(attemptNumber: Int) = CopyingOutputsFailedResponse(helper.jobDescriptorKey, attemptNumber, CopyAttemptError(copyFailureReason))
  def cacheHitBlacklistedResponse(attemptNumber: Int) = CopyingOutputsFailedResponse(helper.jobDescriptorKey, attemptNumber, BlacklistSkip(MetricableCacheCopyErrorCategory.BucketBlacklisted))
}
