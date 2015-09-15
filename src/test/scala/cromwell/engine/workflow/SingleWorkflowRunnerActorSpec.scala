package cromwell.engine.workflow

import cromwell.util.SampleWdl.ThreeStep
import cromwell.{CromwellSpec, CromwellTestkitSpec}

import scala.language.postfixOps

class SingleWorkflowRunnerActorSpec extends CromwellTestkitSpec("SingleWorkflowRunnerActorSpec") {
  val workflowManagerActor = system.actorOf(WorkflowManagerActor.props(dataAccess, CromwellSpec.BackendInstance))
  val props = SingleWorkflowRunnerActor.props(ThreeStep.asWorkflowSources(), ThreeStep.rawInputs, workflowManagerActor)

  "A SingleWorkflowRunnerActor" should {
    "successfully run a workflow" in {
      waitForPattern("workflow finished with status 'Succeeded'.") {
        system.actorOf(props)
      }
    }
  }
}
