package cromwell.engine

import akka.testkit.EventFilter
import cromwell.CromwellTestkitSpec
import cromwell.util.SampleWdl.ThreeStep

import scala.language.postfixOps

class SingleWorkflowRunnerActorSpec extends CromwellTestkitSpec("SingleWorkflowRunnerActorSpec") {
  val workflowManagerActor = this.system.actorOf(WorkflowManagerActor.props)
  val props = SingleWorkflowRunnerActor.props(ThreeStep.WdlSource, ThreeStep.RawInputs, workflowManagerActor)

  "A SingleWorkflowRunnerActor" should {
    "successfully run a workflow" in {
      EventFilter.info(message = "SingleWorkflowRunnerActor: workflow finished with status 'Succeeded'", occurrences = 1) intercept {
        system.actorOf(props)
      }
    }
  }
}
