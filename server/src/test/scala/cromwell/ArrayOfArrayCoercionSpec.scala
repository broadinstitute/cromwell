package cromwell

import akka.testkit.EventFilter
import cromwell.util.SampleWdl
import wom.types._
import wom.values._

class ArrayOfArrayCoercionSpec extends CromwellTestKitWordSpec {
  "A workflow that has an Array[Array[File]] input " should {
    "accept an Array[Array[String]] as the value for the input" in {
      runWdlAndAssertOutputs(
        sampleWdl = SampleWdl.ArrayOfArrays,
        eventFilter = EventFilter.info(pattern = "Workflow complete", occurrences = 1),
        expectedOutputs = Map(
          "wf.subtask.concatenated" -> WomArray(WomArrayType(WomStringType), Seq(
            WomString("foo\nbar\nbaz"),
            WomString("third\nfourth")
          ))
        )
      )
    }
  }
}
