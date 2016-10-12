package cromwell.engine

import cromwell.CromwellTestKitSpec
import cromwell.engine.workflow.WorkflowDescriptorBuilder
import cromwell.util.SampleWdl


class WorkflowManagerActorSpec extends CromwellTestKitSpec with WorkflowDescriptorBuilder {
  override implicit val actorSystem = system

  "A WorkflowManagerActor" should {

    "run workflows in the correct directory" in {
      val outputs = runWdl(sampleWdl = SampleWdl.CurrentDirectory)

      val outputName = "wf_whereami_whereami_pwd"
      val salutation = outputs(outputName)
      val actualOutput = salutation.valueString.trim
      actualOutput should endWith("/call-whereami/execution")
    }
  }
}
