package cromwell.engine.workflow

import akka.testkit.TestActorRef
import cromwell.CromwellTestkitSpec
import cromwell.engine.backend.local.LocalBackend
import cromwell.util.SampleWdl.ThreeStep

import scala.language.postfixOps

class SingleWorkflowRunnerActorSpec extends CromwellTestkitSpec("SingleWorkflowRunnerActorSpec") {
  val workflowManager = TestActorRef(new WorkflowManagerActor(new LocalBackend))
  val props = SingleWorkflowRunnerActor.props(ThreeStep.asWorkflowSources(), ThreeStep.rawInputs, workflowManager)

  "A SingleWorkflowRunnerActor" should {
    "successfully run a workflow" in {
      waitForPattern("workflow finished with status 'Succeeded'.") {
        system.actorOf(props)
      }
    }
  }
}
