package cromwell.engine.workflow.lifecycle.execution.ejea

import cromwell.backend.BackendJobExecutionActor.{AbortedResponse, JobFailedNonRetryableResponse, JobFailedRetryableResponse, JobSucceededResponse}
import cromwell.core.JobOutput
import cromwell.core.callcaching._
import cromwell.engine.workflow.lifecycle.execution.EngineJobExecutionActor.{EJEAData, SucceededResponseData, UpdatingCallCache, UpdatingJobStore}
import cromwell.engine.workflow.lifecycle.execution.callcaching.EngineJobHashingActor.CallCacheHashes
import cromwell.engine.workflow.lifecycle.execution.callcaching.CallCachingEntryId
import cromwell.jobstore.JobStoreActor.RegisterJobCompleted
import cromwell.jobstore.{JobResultSuccess, JobStoreKey}
import org.scalatest.concurrent.Eventually
import wdl4s.values.{WdlInteger, WdlString}

import scala.util.Success

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
    helper.callCacheWriteActorCreations match {
      case GotOne(creation) =>
        creation._1 should be(expectedCallCacheHashes)
        creation._2 should be(expectedResponse)
      case _ => fail("Expected exactly one cache write actor creation.")
    }
    ()
  }
}

private[ejea] trait CanExpectJobStoreWrites extends CanValidateJobStoreKey { self: EngineJobExecutionActorSpec with HasJobSuccessResponse =>

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
}

private[ejea] trait CanExpectHashingInitialization extends Eventually { self: EngineJobExecutionActorSpec =>
  def expectHashingActorInitialization(mode: CallCachingMode): Unit = {
    eventually { helper.jobHashingInitializations.hasExactlyOne should be(true) }
    helper.jobHashingInitializations.checkIt { initialization =>
      initialization._1 should be(helper.backendJobDescriptor)
      initialization._2 should be(mode)
    }
  }
}

private[ejea] trait CanExpectFetchCachedResults extends Eventually { self: EngineJobExecutionActorSpec =>
  def expectFetchCachedResultsActor(expectedCallCachingEntryId: CallCachingEntryId): Unit = {
    eventually { helper.fetchCachedResultsActorCreations.hasExactlyOne should be(true) }
    helper.fetchCachedResultsActorCreations.checkIt {
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
  val successOutputs = Map("a" -> JobOutput(WdlInteger(3)), "b" -> JobOutput(WdlString("bee")))
  def successResponse = JobSucceededResponse(helper.jobDescriptorKey, successRc, successOutputs, None, Seq.empty)
}
private[ejea] object HasJobSuccessResponse {
  val SuccessfulCallCacheHashes = CallCacheHashes(Set(HashResult(HashKey("whatever you want"), HashValue("whatever you need"))))
}

private[ejea] trait HasJobFailureResponses { self: EngineJobExecutionActorSpec =>
  val failedRc = Option(12)
  val failureReason = new Exception("The sixth sheik's sheep is sick!")
  // Need to delay making the response because job descriptors come from the per-test "helper", which is null outside tests!
  def failureRetryableResponse = JobFailedRetryableResponse(helper.jobDescriptorKey, failureReason, failedRc)
  def failureNonRetryableResponse = JobFailedNonRetryableResponse(helper.jobDescriptorKey, failureReason, Option(12))
  def abortedResponse = AbortedResponse(helper.jobDescriptorKey)
}