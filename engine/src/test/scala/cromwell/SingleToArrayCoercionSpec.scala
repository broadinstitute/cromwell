package cromwell

import akka.testkit.EventFilter
import wdl4s.values.{WdlString, WdlFile}
import cromwell.util.SampleWdl

class SingleToArrayCoercionSpec extends CromwellTestkitSpec {
  "A workflow which sends a single output into an Array input" should {
    "be fine, basically" in {
      runWdlAndAssertOutputs(
        sampleWdl = SampleWdl.SingleToArrayCoercion,
        EventFilter.info(pattern = s"Workflow complete", occurrences = 1),
        expectedOutputs = Map(
          "oneToMany.singleFile.out" -> WdlFile("stdout"),
          "oneToMany.listFiles.result" -> WdlString("hello")
        )
      )
    }
  }
}
