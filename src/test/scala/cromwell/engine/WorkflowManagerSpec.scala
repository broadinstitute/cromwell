package cromwell.engine

import java.util.UUID

import cromwell.binding._
import cromwell.engine.backend.local.LocalBackend
import org.scalatest.{Matchers, FlatSpec}
import cromwell.HelloWorldActorSpec._

import scala.util.Try

object TestWorkflowManager {
  val WorkflowString = "hi"
}

class TestWorkflowManager extends WorkflowManager {
  override type Workflow = String
  override val backend = new LocalBackend

  override def generateWorkflow(wdl: WdlSource, inputs: WorkflowInputs): Try[Workflow] = { Try(TestWorkflowManager.WorkflowString)}
  override def workflowStatus(id: WorkflowId): Option[WorkflowState] = {Option(WorkflowRunning)}
}

class WorkflowManagerSpec extends FlatSpec with Matchers {
  "A WorkflowManager" should "allow you to submit a workflow" in {
    val workflowManager = new TestWorkflowManager
    val wf = workflowManager.submitWorkflow(HelloWdl, HelloInputs).get
    wf.id shouldBe a [UUID]
    wf.workflow shouldEqual TestWorkflowManager.WorkflowString
  }

  it should "look up a workflow by id" in {
    val workflowManager = new TestWorkflowManager
    val wf = workflowManager.submitWorkflow(HelloWdl, HelloInputs).get
    workflowManager.workflowById(wf.id).get shouldEqual wf.workflow
  }
}
