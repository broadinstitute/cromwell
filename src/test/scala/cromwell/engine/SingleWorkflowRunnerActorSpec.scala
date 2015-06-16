package cromwell.engine

import cromwell.HelloWorldActorSpec._
import cromwell.util.SampleWdl.ThreeStep
import akka.actor.ActorSystem
import akka.testkit.EventFilter
import com.typesafe.config.ConfigFactory
import cromwell.CromwellTestkitSpec
import scala.language.postfixOps

class SingleWorkflowRunnerActorSpec extends CromwellTestkitSpec(ActorSystem("ActorWorkflowManagerSpec", ConfigFactory.parseString(Config))) {
  val actorSystem = super.getActorSystem
  val workflowManagerActor = actorSystem.actorOf(WorkflowManagerActor.props)
  val props = SingleWorkflowRunnerActor.props(ThreeStep.WdlSource, ThreeStep.RawInputs, workflowManagerActor)

  "A SingleWorkflowRunnerActor" should {
    "successfully run a workflow" in {
      EventFilter.info(message = "Workflow complete: Succeeded", occurrences = 1) intercept {
        actorSystem.actorOf(props)
      }
    }
  }
}