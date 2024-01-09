package cromwell.engine.workflow.lifecycle.execution.ejea

import EngineJobExecutionActorSpec._
import cromwell.backend.BackendJobExecutionActor._
import cromwell.engine.workflow.lifecycle.execution.job.EngineJobExecutionActor._
import cromwell.jobstore.JobStoreActor.{JobStoreWriteFailure, JobStoreWriteSuccess}
import cromwell.engine.workflow.lifecycle.execution.ejea.HasJobSuccessResponse.SuccessfulCallCacheHashes

import scala.util.Success

class EjeaUpdatingJobStoreSpec
    extends EngineJobExecutionActorSpec
    with HasJobSuccessResponse
    with HasJobFailureResponses {

  implicit override val stateUnderTest = UpdatingJobStore

  "An EJEA in UpdatingJobStoreSpec" should {

    List(
      ("SucceededResponse", () => successResponse, true),
      ("FailedRetryableResponse", () => failureRetryableResponse, true),
      ("FailedNonRetryableResponse", () => failureNonRetryableResponse, false)
    ) foreach { case (name, responseMaker, retryable @ _) =>
      s"Forward a saved $name response on and shut down, if the JobStore write is successful" in {
        val response = responseMaker.apply()
        ejea = ejeaInUpdatingJobStoreState(response)
        ejea ! JobStoreWriteSuccess(null) // This value's not read...
        helper.replyToProbe.expectMsg(awaitTimeout, response)
        helper.deathwatch.expectTerminated(ejea, awaitTimeout)
      }
    }

    s"Create a suitable failure if the JobStore write fails" in {
      val response = successResponse
      ejea = ejeaInUpdatingJobStoreState(response)
      val exception = new Exception(
        "I loved Ophelia: forty thousand brothers\\ Could not, with all their quantity of love,\\ Make up my sum. What wilt thou do for her?"
      )
      ejea ! JobStoreWriteFailure(exception)
      helper.replyToProbe.expectMsgPF(awaitTimeout) {
        case JobFailedNonRetryableResponse(jobDescriptorKey, reason, None) =>
          jobDescriptorKey should be(helper.jobDescriptorKey)
          reason.getCause should be(exception)
      }
      helper.deathwatch.expectTerminated(ejea, awaitTimeout)
    }
  }

  def ejeaInUpdatingJobStoreState(response: BackendJobExecutionResponse) = {
    val pendingResponseData =
      ResponsePendingData(helper.backendJobDescriptor, helper.bjeaProps, Some(Success(SuccessfulCallCacheHashes)))
    val newData = response match {
      case success: JobSucceededResponse => pendingResponseData.withSuccessResponse(success)
      case failed: BackendJobFailedResponse => pendingResponseData.withFailedResponse(failed)
      case aborted: JobAbortedResponse => pendingResponseData.withAbortedResponse(aborted)
    }
    helper.buildEJEA().setStateInline(data = newData)
  }

}
