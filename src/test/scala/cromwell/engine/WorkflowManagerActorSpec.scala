package cromwell.engine

import java.util.{Calendar, UUID}

import akka.testkit.TestActorRef
import cromwell.binding.values.WdlString
import cromwell.engine.db.{QueryWorkflowExecutionResult, DummyDataAccess}
import cromwell.engine.workflow.WorkflowManagerActor
import cromwell.engine.workflow.WorkflowManagerActor.{SubmitWorkflow, WorkflowOutputs, WorkflowStatus}
import cromwell.util.SampleWdl
import cromwell.util.SampleWdl.HelloWorld
import cromwell.{CromwellTestkitSpec, binding}


class WorkflowManagerActorSpec extends CromwellTestkitSpec("WorkflowManagerActorSpec") {

  "An ActorWorkflowManager" should {
    "run the Hello World workflow" in {
      implicit val workflowManagerActor = TestActorRef(WorkflowManagerActor.props(DummyDataAccess()), self, "Test the WorkflowManagerActor")

      val workflowId = waitForHandledMessagePattern(pattern = "Transition\\(.*,Running,Succeeded\\)$") {
        messageAndWait[WorkflowId](SubmitWorkflow(HelloWorld.WdlSource, HelloWorld.RawInputs))
      }

      val status = messageAndWait[Option[WorkflowState]](WorkflowStatus(workflowId)).get
      status shouldEqual WorkflowSucceeded

      val outputs = messageAndWait[binding.WorkflowOutputs](WorkflowOutputs(workflowId))

      val actual = outputs.map { case (k, WdlString(string)) => k -> string }
      actual shouldEqual Map(HelloWorld.OutputKey -> HelloWorld.OutputValue)
    }
  }

  it should {
    "Not try to restart any workflows when there are no workflows in restartable states" in {
      waitForPattern("Found no workflows to restart.") {
        TestActorRef(WorkflowManagerActor.props(DummyDataAccess()), self, "No workflows")
      }
    }
  }

  def result(workflowState: WorkflowState,
             wdlSource: String = SampleWdl.HelloWorld.WdlSource,
             wdlInputs: String = SampleWdl.HelloWorld.JsonInputs): QueryWorkflowExecutionResult = {
    QueryWorkflowExecutionResult(
      UUID.randomUUID(), "http://wdl.me", workflowState, Calendar.getInstance().getTime, None, Set.empty, Set.empty, wdlSource, wdlInputs)
  }

  "An WorkflowManagerActor" should {
    "Try to restart workflows when there are workflows in restartable states" in {
      val workflows = Seq(result(WorkflowSubmitted), result(WorkflowRunning))
      val ids = workflows.map { _.workflowId.toString }.sorted
      val dataAccess = new DummyDataAccess {
        override def queryHelper() = workflows
      }
      waitForPattern("Starting workflow IDs: "  + ids.mkString(", ")) {
        waitForPattern("Found 2 workflows to restart.") {
          // TODO this isn't mocking out any part of the restart logic of the WorkflowManagerActor,
          // TODO it actually does "restart" these two hello world workflows.
          TestActorRef(WorkflowManagerActor.props(dataAccess), self, "2 restartable workflows")
        }
      }
    }
  }
}
