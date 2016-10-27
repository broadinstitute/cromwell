package cromwell.engine

import cromwell.CromwellTestkitSpec
import cromwell.engine.workflow.WorkflowDescriptorBuilder
import cromwell.util.SampleWdl


class WorkflowManagerActorSpec extends CromwellTestkitSpec with WorkflowDescriptorBuilder {
  override implicit val actorSystem = system

  "A WorkflowManagerActor" should {

    "run workflows in the correct directory" in {
      val outputs = runWdl(sampleWdl = SampleWdl.CurrentDirectory)

      val outputName = "wf_whereami.whereami.pwd"
      val salutation = outputs(outputName)
      val actualOutput = salutation.valueString.trim
      actualOutput should endWith("/call-whereami/execution")
    }
  }
}
