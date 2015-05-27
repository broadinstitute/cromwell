package cromwell.engine

import java.util.UUID

import cromwell.binding
import cromwell.binding._
import cromwell.engine.backend.local.LocalBackend
import org.scalatest.concurrent.ScalaFutures
import org.scalatest.{Matchers, FlatSpec}
import cromwell.HelloWorldActorSpec._
import scala.concurrent.ExecutionContext.Implicits.global
import scala.concurrent.Future
import scala.util.{Success, Try}

object TestWorkflowManager {
  val WorkflowString = "hi"
}

class TestWorkflowManager extends WorkflowManager {
  override type Workflow = String
  override val backend = new LocalBackend

  override def generateWorkflow(id: WorkflowId, wdl: WdlSource, inputs: WorkflowRawInputs): Try[Workflow] = Try(TestWorkflowManager.WorkflowString)
  override def workflowStatus(id: WorkflowId): Option[WorkflowState] = Option(WorkflowRunning) // FIXME
  override def workflowOutputs(id: WorkflowId): Future[Option[binding.WorkflowOutputs]] = {Future {Option(Map.empty)}} // FIXME
  override def submitWorkflow(wdl: WdlSource, inputs: binding.WorkflowRawInputs): Future[Try[WorkflowId]] = Future {Success(UUID.randomUUID())}
  override def updateWorkflowState(workflow: Workflow, state: WorkflowState): Unit = {} // FIXME: This doesn't do anything
}

class WorkflowManagerSpec extends FlatSpec with Matchers with ScalaFutures {
  "A WorkflowManager" should "allow you to submit a workflow" in {
    val workflowManager = new TestWorkflowManager
    val eventualSubmission = workflowManager.submitWorkflow(HelloWdl, HelloRawInputs)
    whenReady(eventualSubmission) { w =>
      val wfId = w.get
      wfId shouldBe a [UUID]
    }
  }
}
