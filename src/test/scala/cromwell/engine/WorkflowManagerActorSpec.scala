package cromwell.engine

import java.util.{Calendar, UUID}

import akka.testkit.TestActorRef
import cromwell.binding.FullyQualifiedName
import cromwell.binding.types.WdlStringType
import cromwell.binding.values.WdlString
import cromwell.engine.ExecutionStatus.{Done, NotStarted, Running}
import cromwell.engine.db.DataAccess.WorkflowInfo
import cromwell.engine.db.{CallStatus, QueryWorkflowExecutionResult, DummyDataAccess}
import cromwell.engine.workflow.WorkflowManagerActor
import cromwell.engine.workflow.WorkflowManagerActor.{SubmitWorkflow, WorkflowOutputs, WorkflowStatus}
import cromwell.util.SampleWdl
import cromwell.util.SampleWdl.HelloWorld
import cromwell.{engine, CromwellTestkitSpec, binding}

import scala.concurrent.Future


class WorkflowManagerActorSpec extends CromwellTestkitSpec("WorkflowManagerActorSpec") {

  "A WorkflowManagerActor" should {

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
      val ids = workflows.map {
        _.workflowId.toString
      }.sorted
      val key = SymbolStoreKey("hello.hello", "addressee", None, input = true)
      val symbols = Map(key -> new SymbolStoreEntry(key, WdlStringType, Option(WdlString("world"))))

      val dataAccess = new DummyDataAccess() {
        workflows foreach { workflow =>
          val id = workflow.workflowId
          executionStatuses += (id -> Map("hello.hello" -> (if (id == submitted.workflowId) NotStarted else Running)))
          symbolStore += (workflow.workflowId -> symbols)
        }

        override def getWorkflowsByState(states: Traversable[WorkflowState]): Future[Traversable[WorkflowInfo]] = {
          Future.successful {
            workflows.map { w =>
              WorkflowInfo(w.workflowId, w.wdlSource, w.jsonInputs)
            }
          }
        }
      }

      waitForPattern("Restarting workflow IDs: " + ids.mkString(", ")) {
        waitForPattern("Found 2 workflows to restart.") {
          // Workflows are always set back to Submitted on restart.
          waitForPattern("transitioning from Submitted to Running.", occurrences = 2) {
            // Both the previously in-flight call and the never-started call should get started.
            waitForPattern("Starting calls: hello.hello", occurrences = 2) {
              waitForPattern("transitioning from Running to Succeeded", occurrences = 2) {
                TestActorRef(WorkflowManagerActor.props(dataAccess), self, "2 restartable workflows")
              }
            }
          }
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
