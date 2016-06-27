package cromwell

import akka.testkit.{EventFilter, TestActorRef}
import cromwell.backend.PreemptedException
import cromwell.core.WorkflowSucceeded
import cromwell.engine._
import cromwell.engine.backend._
import cromwell.engine.workflow.WorkflowDescriptorBuilder
import cromwell.util.SampleWdl
import cromwell.webservice.CromwellApiHandler._
import org.specs2.mock.Mockito
import wdl4s.values.{WdlArray, WdlValue}

import scala.concurrent.{ExecutionContext, Future}
import scala.language.postfixOps


class RetryableCallsSpec extends CromwellTestkitSpec with Mockito with WorkflowDescriptorBuilder {
  override val actorSystem = system
//  val customizedLocalBackend = new OldStyleLocalBackend(CromwellTestkitSpec.DefaultLocalBackendConfigEntry, system) {
//    override def execute(jobDescriptor: OldStyleBackendCallJobDescriptor)(implicit ec: ExecutionContext): Future[OldStyleExecutionHandle] = {
//      jobDescriptor.key.scope.taskFqn match {
//        case "do_scatter" if jobDescriptor.key.index.contains(0) && jobDescriptor.key.attempt == 1 =>
//          RetryableExecutionHandle(new PreemptedException("Retryable failure")).future
//        case _ => super.execute(jobDescriptor)
//      }
//    }
//  }

  val customizedLocalBackendOptions =
    """
      |{
      |  "backend": "customizedLocalBackend"
      |}""".stripMargin

  "A workflow with a Successful retried Call" should {

    "succeed and return correct metadata" ignore {
//      CromwellBackends.registerCustomBackend("customizedLocalBackend", customizedLocalBackend)
//      implicit val workflowManagerActor = TestActorRef(OldStyleWorkflowManagerActor.props(), self, "Test Workflow metadata with Retried Calls")
//      val wd = materializeWorkflowDescriptorFromSources(workflowSources = SampleWdl.PrepareScatterGatherWdl().asWorkflowSources(workflowOptions = customizedLocalBackendOptions))
//      val workflowId = waitForHandledMessagePattern(pattern = "transitioning from Running to Succeeded") {
//        EventFilter.info(pattern = s"persisting status of do_scatter:0:2 to Starting", occurrences = 1).intercept {
//          messageAndWait[WorkflowManagerSubmitSuccess](SubmitWorkflow(SampleWdl.PrepareScatterGatherWdl().asWorkflowSources(workflowOptions = customizedLocalBackendOptions))).id
//        }
//      }
//
//      val status = messageAndWait[WorkflowManagerStatusSuccess](WorkflowStatus(workflowId)).state
//      status shouldEqual WorkflowSucceeded
//
//      val metadata = messageAndWait[WorkflowManagerWorkflowMetadataSuccess](WorkflowMetadata(workflowId)).response
//      metadata should not be null
//
//      metadata.status shouldBe WorkflowSucceeded.toString
//      metadata.start shouldBe defined
//      metadata.end shouldBe defined
//      metadata.outputs shouldBe defined
//      metadata.outputs.get should have size 3
//      metadata.calls should have size 3
//      metadata.calls.get("sc_test.do_scatter") shouldBe defined
//      metadata.calls("sc_test.do_scatter") should have size 5
//
//      metadata.calls("sc_test.do_scatter").head.executionStatus shouldBe "Preempted"
//      metadata.calls("sc_test.do_scatter").head.shardIndex shouldBe 0
//      metadata.calls("sc_test.do_scatter").head.attempt shouldBe 1
//      metadata.calls("sc_test.do_scatter").head.stdout shouldBe defined
//      metadata.calls("sc_test.do_scatter").head.stdout.get.value should not include "attempt"
//      metadata.calls("sc_test.do_scatter").head.stderr shouldBe defined
//      metadata.calls("sc_test.do_scatter").head.stderr.get.value should not include "attempt"
//      metadata.calls("sc_test.do_scatter").head.outputs.get should have size 0
//      metadata.calls("sc_test.do_scatter")(1).executionStatus shouldBe "Done"
//      metadata.calls("sc_test.do_scatter")(1).shardIndex shouldBe 0
//      metadata.calls("sc_test.do_scatter")(1).attempt shouldBe 2
//      metadata.calls("sc_test.do_scatter")(1).stdout shouldBe defined
//      metadata.calls("sc_test.do_scatter")(1).stdout.get.value should include ("attempt")
//      metadata.calls("sc_test.do_scatter")(1).stderr shouldBe defined
//      metadata.calls("sc_test.do_scatter")(1).stderr.get.value should include ("attempt")
//      metadata.calls("sc_test.do_scatter")(1).outputs shouldBe defined
//      metadata.calls("sc_test.do_scatter")(1).outputs.get should not be empty
//
//      // Check that only successful shards (and not preempted ones) are collected
//      val gatherInputs: Seq[WdlValue] = metadata.calls("sc_test.do_gather").head.inputs("input_files").asInstanceOf[WdlArray].value
//      gatherInputs.size shouldBe 4
//
//      val wfStdouterr = messageAndWait[WorkflowManagerWorkflowStdoutStderrSuccess](WorkflowStdoutStderr(workflowId))
//      wfStdouterr should not be null
//      wfStdouterr.logs should have size 3
//      wfStdouterr.logs("sc_test.do_scatter") should have size 5
//      wfStdouterr.logs("sc_test.do_scatter").head.stdout.value should not include "attempt"
//      wfStdouterr.logs("sc_test.do_scatter").head.stderr.value should not include "attempt"
//      wfStdouterr.logs("sc_test.do_scatter")(1).stdout.value should include ("attempt")
//      wfStdouterr.logs("sc_test.do_scatter")(1).stderr.value should include ("attempt")
//
//      val callStdouterr = messageAndWait[WorkflowManagerCallStdoutStderrSuccess](CallStdoutStderr(workflowId, "sc_test.do_scatter"))
//      callStdouterr should not be null
//      callStdouterr.logs should have size 5
//      callStdouterr.logs.head.stdout.value should not include "attempt"
//      callStdouterr.logs.head.stderr.value should not include "attempt"
//      callStdouterr.logs(1).stdout.value should include ("attempt")
//      callStdouterr.logs(1).stderr.value should include ("attempt")

    }

  }
}
