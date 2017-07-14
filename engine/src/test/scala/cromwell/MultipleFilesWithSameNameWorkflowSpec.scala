package cromwell

import akka.testkit._
import cromwell.util.SampleWdl
import wdl4s.wdl.values.WdlString


class MultipleFilesWithSameNameWorkflowSpec extends CromwellTestKitWordSpec {
  "A workflow with two file inputs that have the same name" should {
    "not clobber one file with the contents of another" in {
      runWdlAndAssertOutputs(
        sampleWdl = SampleWdl.FileClobber,
        EventFilter.info(pattern = "Starting calls: two.x:NA:1, two.y:NA:1", occurrences = 1),
        expectedOutputs = Map(
          "two.x.out" -> WdlString("first file.txt"),
          "two.y.out" -> WdlString("second file.txt")
        )
      )
    }
  }
}
