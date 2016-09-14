package cromwell.engine.workflow.lifecycle.execution.ejea

import akka.testkit.TestProbe
import cromwell.backend.BackendJobDescriptorKey
import cromwell.backend.BackendJobExecutionActor.{FailedNonRetryableResponse, FailedRetryableResponse, RecoverJobCommand, SucceededResponse}
import cromwell.core._
import cromwell.engine.workflow.lifecycle.execution.EngineJobExecutionActor.{CheckingJobStore, JobRunning, NoData, PreparingJob}
import cromwell.engine.workflow.lifecycle.execution.JobPreparationActor.BackendJobPreparationFailed
import cromwell.jobstore.{JobResultFailure, JobResultSuccess}
import cromwell.jobstore.JobStoreActor.{JobComplete, JobNotComplete}
import EngineJobExecutionActorSpec.EnhancedTestEJEA
import cromwell.engine.workflow.lifecycle.execution.JobPreparationActor

class EjeaCheckingJobStoreSpec extends EngineJobExecutionActorSpec {

  override implicit val stateUnderTest = CheckingJobStore

  "An EJEA in CheckingJobStore state should" should {
    "send a Job SucceededResponse if the job is already complete and successful" in {
      createCheckingJobStoreEjea()
      ejea.setState(CheckingJobStore)
      val returnCode: Option[Int] = Option(0)
      val jobOutputs: JobOutputs = Map.empty

      ejea ! JobComplete(JobResultSuccess(returnCode, jobOutputs))

      helper.replyToProbe.expectMsgPF(awaitTimeout) {
        case response: SucceededResponse =>
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
          case response: FailedNonRetryableResponse =>
            false should be(retryable)
            response.returnCode shouldBe returnCode
            response.throwable shouldBe reason
          case response: FailedRetryableResponse =>
            true should be(retryable)
            response.returnCode shouldBe returnCode
            response.throwable shouldBe reason
        }

        helper.deathwatch.expectTerminated(ejea)
      }
    }

    "begin preparing the job if it's not already complete" in {
      createCheckingJobStoreEjea()
      ejea.setState(CheckingJobStore)
      ejea ! JobNotComplete

      helper.jobPreparationProbe.expectMsg(awaitTimeout, "expecting RecoverJobCommand", JobPreparationActor.Start)
      ejea.stateName should be(PreparingJob)

      ejea.stop()
    }
  }

  private def createCheckingJobStoreEjea(): Unit = { ejea = helper.buildEJEA(restarting = true).setStateInline(state = CheckingJobStore, data = NoData) }
}
