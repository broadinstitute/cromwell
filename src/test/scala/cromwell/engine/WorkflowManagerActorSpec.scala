package cromwell.engine

import java.util.{Calendar, UUID}

import akka.testkit.TestActorRef
import cromwell.binding.FullyQualifiedName
import cromwell.binding.values.WdlString
import cromwell.engine.db.DataAccess.WorkflowInfo
import cromwell.engine.db.{CallStatus, QueryWorkflowExecutionResult, DummyDataAccess}
import cromwell.engine.workflow.WorkflowManagerActor
import cromwell.engine.workflow.WorkflowManagerActor.{SubmitWorkflow, WorkflowOutputs, WorkflowStatus}
import cromwell.util.SampleWdl
import cromwell.util.SampleWdl.HelloWorld
import cromwell.{CromwellTestkitSpec, binding}

import scala.concurrent.Future


class WorkflowManagerActorSpec extends CromwellTestkitSpec("WorkflowManagerActorSpec") {

  "An ActorWorkflowManager" should {

    "run the Hello World workflow" in {
      implicit val workflowManagerActor = TestActorRef(WorkflowManagerActor.props(DummyDataAccess()), self, "Test the WorkflowManagerActor")

      val workflowId = waitForHandledMessagePattern(pattern = "Transition\\(.*,Running,Succeeded\\)$") {
        messageAndWait[WorkflowId](SubmitWorkflow(HelloWorld.WdlSource, HelloWorld.WdlJson, HelloWorld.RawInputs))
      }

      val status = messageAndWait[Option[WorkflowState]](WorkflowStatus(workflowId)).get
      status shouldEqual WorkflowSucceeded

      val outputs = messageAndWait[binding.WorkflowOutputs](WorkflowOutputs(workflowId))

      val actual = outputs.map { case (k, WdlString(string)) => k -> string }
      actual shouldEqual Map(HelloWorld.OutputKey -> HelloWorld.OutputValue)
    }

    "Not try to restart any workflows when there are no workflows in restartable states" in {
      waitForPattern("Found no workflows to restart.") {
        TestActorRef(WorkflowManagerActor.props(DummyDataAccess()), self, "No workflows")
      }
    }

    "Try to restart workflows when there are workflows in restartable states" in {
      val (submitted, running) = (result(WorkflowSubmitted), result(WorkflowRunning))
      val workflows = Seq(submitted, running)
      val ids = workflows.map { _.workflowId.toString }.sorted
      val dataAccess = new DummyDataAccess() {

        override def getWorkflowsByState(states: Traversable[WorkflowState]): Future[Traversable[WorkflowInfo]] = {
          Future.successful { workflows.map { w =>
              WorkflowInfo(w.workflowId, w.wdlSource, w.jsonInputs)
            }
          }
        }

        override def getExecutionStatuses(workflowId: WorkflowId): Future[Map[FullyQualifiedName, CallStatus]] = {
          Future.successful { workflows.map { w =>
            "hello" -> (if (w.workflowId == submitted.workflowId) ExecutionStatus.NotStarted else ExecutionStatus.Running) }.toMap
          }
        }
      }
      waitForPattern("Restarting workflow IDs: "  + ids.mkString(", ")) {
        waitForPattern("Found 2 workflows to restart.") {
          // TODO this isn't mocking out any part of the restart logic of the WorkflowManagerActor,
          // TODO it actually does "restart" these two hello world workflows.
          TestActorRef(WorkflowManagerActor.props(dataAccess), self, "2 restartable workflows")
        }
      }
    }
  }

  def result(workflowState: WorkflowState,
             wdlSource: String = SampleWdl.HelloWorld.WdlSource,
             wdlInputs: String = SampleWdl.HelloWorld.WdlJson): QueryWorkflowExecutionResult = {
    QueryWorkflowExecutionResult(
      UUID.randomUUID(), "http://wdl.me", workflowState, Calendar.getInstance().getTime, None, Set.empty, Set.empty, wdlSource, wdlInputs)
  }
}
