package cromwell.engine

import akka.testkit._
import cromwell.CromwellSpec.PostMVP
import cromwell.CromwellTestkitSpec
import cromwell.engine.workflow.WorkflowDescriptorBuilder
import cromwell.util.SampleWdl

import scala.concurrent.duration._
import scala.language.postfixOps

class WorkflowManagerActorSpec extends CromwellTestkitSpec with WorkflowDescriptorBuilder {
  override implicit val actorSystem = system

  "A WorkflowManagerActor" should {

    val TestExecutionTimeout = 5.seconds.dilated

    // TODO PBE: Restart workflows tests: re-add (but somewhere else?) in 0.21
    "Not try to restart any workflows when there are no workflows in restartable states" taggedAs PostMVP ignore {

    }


    "Try to restart workflows when there are workflows in restartable states" taggedAs PostMVP ignore {

    }

    "run workflows in the correct directory" in {
      val outputs = runWdl(sampleWdl = SampleWdl.CurrentDirectory)

      val outputName = "whereami.whereami.pwd"
      val salutation = outputs.get(outputName).get
      val actualOutput = salutation.valueString.trim
      actualOutput should endWith("/call-whereami")
    }

  }
}
