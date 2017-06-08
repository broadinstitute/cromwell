package cromwell.engine.workflow.lifecycle.execution.ejea

import cromwell.engine.workflow.lifecycle.execution.EngineJobExecutionActor.{FailedNonRetryableResponseData, ResponsePendingData, RunningJob, SucceededResponseData}
import cromwell.jobstore.JobResultFailure
import cromwell.jobstore.JobStoreActor.RegisterJobCompleted
import EngineJobExecutionActorSpec.EnhancedTestEJEA
import cromwell.core.callcaching._
import cromwell.engine.workflow.lifecycle.execution.callcaching.EngineJobHashingActor.HashError
import org.scalatest.concurrent.Eventually
import cromwell.engine.workflow.lifecycle.execution.ejea.HasJobSuccessResponse.SuccessfulCallCacheHashes

import scala.util.{Failure, Success}

class EjeaRunningJobSpec extends EngineJobExecutionActorSpec with Eventually with CanValidateJobStoreKey with CanExpectJobStoreWrites with CanExpectCacheWrites with HasJobSuccessResponse with HasJobFailureResponses {

  override implicit val stateUnderTest = RunningJob

  val hashError = HashError(new Exception("ARGH!!!"))

  "A 'RunningJob' EJEA" should {
    CallCachingModes foreach { mode =>
      val andMaybeCallCacheHashes = if (mode.writeToCache) "then CallCacheHashes " else ""
      s"Handle receiving a SucceededResponse ${andMaybeCallCacheHashes}correctly in $mode mode" in {
        ejea = ejeaInRunningState(mode)
        ejea ! successResponse

        if (mode.writeToCache) {
          // The data should update. The state shouldn't:
          eventually { ejea.stateData should be(SucceededResponseData(successResponse, None)) }
          ejea.stateName should be(RunningJob)

          // When the hashes arrive, the state and data should update:
          ejea ! SuccessfulCallCacheHashes
          expectCacheWriteForSuccessfulJob(successResponse, SuccessfulCallCacheHashes)
        } else {
          expectJobStoreWrite(SucceededResponseData(successResponse, None))
        }
      }

      if (mode.writeToCache) {
        s"Handle receiving CallCacheHashes then SuccessResponse correctly in $mode mode" in {
          ejea = ejeaInRunningState(mode)
          ejea ! SuccessfulCallCacheHashes
          eventually { ejea.stateData should be(initialData.copy(hashes = Some(Success(SuccessfulCallCacheHashes)))) }
          ejea.stateName should be(RunningJob)
          ejea ! successResponse
          expectCacheWriteForSuccessfulJob(successResponse, SuccessfulCallCacheHashes)
        }

        s"Handle receiving SuccessResponse then HashError correctly in $mode mode" in {
          ejea = ejeaInRunningState(mode)
          ejea ! successResponse
          eventually { ejea.stateData should be(SucceededResponseData(successResponse, None)) }
          ejea.stateName should be(RunningJob)
          ejea ! hashError
          expectJobStoreWrite(SucceededResponseData(successResponse, Some(Failure(hashError.reason))))
        }

        s"Handle receiving HashError then SuccessResponse correctly in $mode mode" in {
          ejea = ejeaInRunningState(mode)
          ejea ! hashError
          eventually { ejea.stateData should be(ResponsePendingData(helper.backendJobDescriptor, helper.bjeaProps, Some(Failure(hashError.reason)))) }
          ejea.stateName should be(RunningJob)
          ejea ! successResponse
          expectJobStoreWrite(SucceededResponseData(successResponse, Some(Failure(hashError.reason))))
        }

        s"Handle receiving CallCacheHashes then FailedNonRetryable correctly in $mode mode" in {
          ejea = ejeaInRunningState(mode)
          ejea ! SuccessfulCallCacheHashes
          eventually { ejea.stateData should be(initialData.copy(hashes = Some(Success(SuccessfulCallCacheHashes)))) }
          ejea.stateName should be(RunningJob)
          ejea ! failureNonRetryableResponse
          expectCacheWriteForFailedNonRetryableJob(failureNonRetryableResponse, SuccessfulCallCacheHashes)
        }

        s"Handle receiving FailedNonRetryable then HashError correctly in $mode mode" in {
          ejea = ejeaInRunningState(mode)
          ejea ! failureNonRetryableResponse
          eventually { ejea.stateData should be(FailedNonRetryableResponseData(failureNonRetryableResponse, None)) }
          ejea.stateName should be(RunningJob)
          ejea ! hashError
          helper.jobStoreProbe.expectMsgPF(max = awaitTimeout, hint = "Job Store Write") {
            case RegisterJobCompleted(jobKey, JobResultFailure(returnCode, reason, isRetryable)) =>
              validateJobStoreKey(jobKey)
              returnCode should be(failedRc)
              reason should be(failureReason)
              isRetryable should be(false)
          }
        }

        s"Handle receiving HashError then FailedNonRetryable correctly in $mode mode" in {
          ejea = ejeaInRunningState(mode)
          ejea ! hashError
          eventually { ejea.stateData should be(ResponsePendingData(helper.backendJobDescriptor, helper.bjeaProps, Some(Failure(hashError.reason)))) }
          ejea.stateName should be(RunningJob)
          ejea ! failureNonRetryableResponse
          helper.jobStoreProbe.expectMsgPF(max = awaitTimeout, hint = "Job Store Write") {
            case RegisterJobCompleted(jobKey, JobResultFailure(returnCode, reason, isRetryable)) =>
              validateJobStoreKey(jobKey)
              returnCode should be(failedRc)
              reason should be(failureReason)
              isRetryable should be(false)
          }
        }
      }
    }

    s"register 'FailedRetryableResponse's with the JobStore" in {
      ejea = ejeaInRunningState()
      val response = failureRetryableResponse
      ejea.underlyingActor.receive.isDefinedAt(response) should be(true)
      ejea ! response

      helper.jobStoreProbe.expectMsgPF(max = awaitTimeout, hint = "Job Store Write") {
        case RegisterJobCompleted(jobKey, JobResultFailure(returnCode, reason, isRetryable)) =>
          validateJobStoreKey(jobKey)
          returnCode should be(failedRc)
          reason should be(failureReason)
          isRetryable should be(true)
      }
    }


    "not register aborted jobs in the job store, forward straight to parent instead" in {
      ejea = ejeaInRunningState()
      ejea.underlyingActor.receive.isDefinedAt(abortedResponse) should be(true)
      ejea ! abortedResponse

      helper.replyToProbe.expectMsg(max = awaitTimeout, hint = "parent wants the response", abortedResponse)

      helper.deathwatch.expectTerminated(ejea, awaitTimeout)
      // Make sure nothing was sent to the JobStore or CacheResultSaver in the meanwhile:
      helper.jobStoreProbe.expectNoMsg(awaitAlmostNothing)
    }
  }

  def initialData = ResponsePendingData(helper.backendJobDescriptor, helper.bjeaProps, None)
  def ejeaInRunningState(mode: CallCachingMode = CallCachingActivity(ReadAndWriteCache)) = helper.buildEJEA(callCachingMode = mode).setStateInline(state = RunningJob, data = initialData)
}
