package cromwell.engine.workflow.lifecycle.execution.ejea

import cromwell.backend.BackendJobExecutionActor.{JobFailedNonRetryableResponse, JobFailedRetryableResponse, JobSucceededResponse}
import cromwell.core._
import cromwell.engine.workflow.lifecycle.execution.EngineJobExecutionActor.{CheckingCacheEntryExistence, CheckingJobStore, NoData}
import cromwell.engine.workflow.lifecycle.execution.callcaching.CallCacheReadActor.CallCacheEntryForCall
import cromwell.engine.workflow.lifecycle.execution.ejea.EngineJobExecutionActorSpec.EnhancedTestEJEA
import cromwell.jobstore.JobStoreActor.{JobComplete, JobNotComplete}
import cromwell.jobstore.{JobResultFailure, JobResultSuccess}

class EjeaCheckingJobStoreSpec extends EngineJobExecutionActorSpec {

  override implicit val stateUnderTest = CheckingJobStore

  "An EJEA in CheckingJobStore state should" should {
    "send a Job SucceededResponse if the job is already complete and successful" in {
      createCheckingJobStoreEjea()
      ejea.setState(CheckingJobStore)
      val returnCode: Option[Int] = Option(0)
      val jobOutputs: CallOutputs = Map.empty

      ejea ! JobComplete(JobResultSuccess(returnCode, jobOutputs))

      helper.replyToProbe.expectMsgPF(awaitTimeout) {
        case response: JobSucceededResponse =>
          response.returnCode shouldBe returnCode
          response.jobOutputs shouldBe jobOutputs
      }

      helper.deathwatch.expectTerminated(ejea)
    }

    List(("FailedNonRetryableResponse", false), ("FailedRetryableResponse", true)) foreach { case (name, retryable) =>

      s"send a $name if the job is already complete and failed" in {
        createCheckingJobStoreEjea()
        val returnCode: Option[Int] = Option(1)
        val reason: Throwable = new Exception("something horrible happened...")

        ejea ! JobComplete(JobResultFailure(returnCode, reason, retryable))

        helper.replyToProbe.expectMsgPF(awaitTimeout) {
          case response: JobFailedNonRetryableResponse =>
            false should be(retryable)
            response.returnCode shouldBe returnCode
            response.throwable shouldBe reason
          case response: JobFailedRetryableResponse =>
            true should be(retryable)
            response.returnCode shouldBe returnCode
            response.throwable shouldBe reason
        }

        helper.deathwatch.expectTerminated(ejea)
      }
    }

    "check for cache entry existence if it's not already complete" in {
      createCheckingJobStoreEjea()
      ejea.setState(CheckingJobStore)
      ejea ! JobNotComplete

      helper.callCacheReadActorProbe.expectMsg(awaitTimeout, "expecting CallCacheEntryForCall", CallCacheEntryForCall(helper.workflowId, helper.jobDescriptorKey))
      ejea.stateName should be(CheckingCacheEntryExistence)

      ejea.stop()
    }
  }

  private def createCheckingJobStoreEjea(): Unit = { ejea = helper.buildEJEA(restarting = true).setStateInline(state = CheckingJobStore, data = NoData) }
}
