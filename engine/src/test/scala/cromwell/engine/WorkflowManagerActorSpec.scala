package cromwell.engine

import cromwell.CromwellTestkitSpec
import cromwell.engine.workflow.WorkflowDescriptorBuilder
import cromwell.util.SampleWdl

import scala.language.postfixOps

class WorkflowManagerActorSpec extends CromwellTestkitSpec with WorkflowDescriptorBuilder {
  override implicit val actorSystem = system

  "A WorkflowManagerActor" should {

    "run workflows in the correct directory" in {
      val outputs = runWdl(sampleWdl = SampleWdl.CurrentDirectory)

      val outputName = "whereami.whereami.pwd"
      val salutation = outputs.get(outputName).get
      val actualOutput = salutation.valueString.trim
      actualOutput should endWith("/call-whereami/execution")
    }
  }
}
