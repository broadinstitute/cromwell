package cromwell.engine.workflow

import akka.testkit.EventFilter
import cromwell.CromwellTestkitSpec
import cromwell.engine.db.DummyDataAccess
import cromwell.util.SampleWdl.ThreeStep

import scala.language.postfixOps

class SingleWorkflowRunnerActorSpec extends CromwellTestkitSpec("ActorWorkflowManagerSpec") {
  val workflowManagerActor = system.actorOf(WorkflowManagerActor.props(DummyDataAccess()))
  val props = SingleWorkflowRunnerActor.props(ThreeStep.WdlSource, ThreeStep.RawInputs, workflowManagerActor)

  "A SingleWorkflowRunnerActor" should {
    "successfully run a workflow" in {
      EventFilter.info(message = "Workflow complete: Succeeded", occurrences = 1) intercept {
        system.actorOf(props)
      }
    }
  }
}