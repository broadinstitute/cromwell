package cromwell.engine.workflow

import cromwell.CromwellTestkitSpec
import cromwell.binding.FullyQualifiedName
import cromwell.engine.ExecutionStatus.NotStarted
import cromwell.engine.WorkflowId
import cromwell.engine.db.{CallStatus, DummyDataAccess}
import cromwell.util.SampleWdl.ThreeStep

import scala.concurrent.Future
import scala.language.postfixOps

class SingleWorkflowRunnerActorSpec extends CromwellTestkitSpec("SingleWorkflowRunnerActorSpec") {
  val dataAccess = new DummyDataAccess() {
    override def getExecutionStatuses(workflowId: WorkflowId): Future[Map[FullyQualifiedName, CallStatus]] = {
      Future.successful(Seq("ps", "cgrep", "wc").map { _ -> NotStarted}.toMap)
    }
  }
  val workflowManagerActor = system.actorOf(WorkflowManagerActor.props(dataAccess))
  val props = SingleWorkflowRunnerActor.props(ThreeStep.WdlSource, ThreeStep.WdlJson, ThreeStep.RawInputs, workflowManagerActor)

  "A SingleWorkflowRunnerActor" should {
    "successfully run a workflow" in {
      waitForPattern("workflow finished with status 'Succeeded'.") {
        system.actorOf(props)
      }
    }
  }
}