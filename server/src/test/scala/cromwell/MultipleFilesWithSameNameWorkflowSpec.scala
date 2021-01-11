package cromwell

import cromwell.util.SampleWdl
import wom.values.WomString

class MultipleFilesWithSameNameWorkflowSpec extends CromwellTestKitWordSpec {
  "A workflow with two file inputs that have the same name" should {
    "not clobber one file with the contents of another" in {
      runWdlAndAssertOutputs(
        sampleWdl = SampleWdl.FileClobber,
        expectedOutputs = Map(
          "two.x.out" -> WomString("first file.txt"),
          "two.y.out" -> WomString("second file.txt")
        ),
        testActorName = "TestCromwellRootActor-clobber",
      )
    }
  }
}
