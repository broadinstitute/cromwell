package cromwell

import akka.testkit.{EventFilter, TestActorRef}
import cromwell.engine.backend._
import cromwell.engine.backend.local.LocalBackend
import cromwell.engine.workflow.WorkflowManagerActor
import cromwell.engine.workflow.WorkflowManagerActor._
import cromwell.engine.{PreemptedException, WorkflowId, WorkflowSucceeded}
import cromwell.util.SampleWdl
import cromwell.webservice.CromwellApiHandler._
import org.specs2.mock.Mockito

import scala.concurrent.{ExecutionContext, Future}

class RetryableCallsSpec extends CromwellTestkitSpec with Mockito {
  val customizedLocalBackend = new LocalBackend(system) {
    var nbCall = 0
    override def execute(backendCall: BackendCall)(implicit ec: ExecutionContext): Future[ExecutionHandle] = {
      backendCall.call.task.name match {
        case "ps" if nbCall < 1 =>
          nbCall += 1
          RetryableExecutionHandle(new PreemptedException("Retryable failure")).future
        case _ => super.execute(backendCall)
      }
    }
  }

  "A workflow with a Successful retried Call" should {

    "succeed and return correct metadata" in {
      implicit val workflowManagerActor = TestActorRef(WorkflowManagerActor.props(customizedLocalBackend), self, "Test Workflow metadata with Retried Calls")

      val workflowId = waitForHandledMessagePattern(pattern = "transitioning from Running to Succeeded") {
        EventFilter.info(pattern = s"starting calls: three_step.ps", occurrences = 2).intercept {
          messageAndWait[WorkflowManagerSubmitSuccess](SubmitWorkflow(SampleWdl.ThreeStep.asWorkflowSources())).id
        }
      }

      val status = messageAndWait[WorkflowManagerStatusSuccess](WorkflowStatus(workflowId)).state
      status shouldEqual WorkflowSucceeded

      val metadata = messageAndWait[WorkflowManagerWorkflowMetadataSuccess](WorkflowMetadata(workflowId)).response
      metadata should not be null

      metadata.status shouldBe WorkflowSucceeded.toString
      metadata.start shouldBe defined
      metadata.end shouldBe defined
      metadata.outputs shouldBe defined
      metadata.outputs.get should have size 3
      metadata.calls should have size 3
      metadata.calls.get("three_step.ps") shouldBe defined
      metadata.calls("three_step.ps") should have size 2

      metadata.calls("three_step.ps").head.executionStatus shouldBe "Preempted"
      metadata.calls("three_step.ps").head.attempt shouldBe 1
      metadata.calls("three_step.ps").head.stdout shouldBe defined
      metadata.calls("three_step.ps").head.stdout.get.value should not include "attempt"
      metadata.calls("three_step.ps").head.stderr shouldBe defined
      metadata.calls("three_step.ps").head.stderr.get.value should not include "attempt"
      metadata.calls("three_step.ps").head.outputs.get should have size 0
      metadata.calls("three_step.ps")(1).executionStatus shouldBe "Done"
      metadata.calls("three_step.ps")(1).attempt shouldBe 2
      metadata.calls("three_step.ps")(1).stdout shouldBe defined
      metadata.calls("three_step.ps")(1).stdout.get.value should include ("attempt")
      metadata.calls("three_step.ps")(1).stderr shouldBe defined
      metadata.calls("three_step.ps")(1).stderr.get.value should include ("attempt")
      metadata.calls("three_step.ps")(1).outputs shouldBe defined
      metadata.calls("three_step.ps")(1).outputs.get should not be empty

      val wfStdouterr = messageAndWait[WorkflowManagerWorkflowStdoutStderrSuccess](WorkflowStdoutStderr(workflowId))
      wfStdouterr should not be null
      wfStdouterr.logs should have size 3
      wfStdouterr.logs("three_step.ps") should have size 2
      wfStdouterr.logs("three_step.ps").head.stdout.value should not include "attempt"
      wfStdouterr.logs("three_step.ps").head.stderr.value should not include "attempt"
      wfStdouterr.logs("three_step.ps")(1).stdout.value should include ("attempt")
      wfStdouterr.logs("three_step.ps")(1).stderr.value should include ("attempt")

      val callStdouterr = messageAndWait[WorkflowManagerCallStdoutStderrSuccess](CallStdoutStderr(workflowId, "three_step.ps"))
      callStdouterr should not be null
      callStdouterr.logs should have size 2
      callStdouterr.logs.head.stdout.value should not include "attempt"
      callStdouterr.logs.head.stderr.value should not include "attempt"
      callStdouterr.logs(1).stdout.value should include ("attempt")
      callStdouterr.logs(1).stderr.value should include ("attempt")

    }

  }
}
