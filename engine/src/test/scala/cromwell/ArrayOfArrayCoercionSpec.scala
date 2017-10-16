package cromwell

import akka.testkit._
import cromwell.util.SampleWdl
import wdl.types.{WdlArrayType, WdlStringType}
import wdl.values.{WdlArray, WdlString}


class ArrayOfArrayCoercionSpec extends CromwellTestKitWordSpec {
  "A workflow that has an Array[Array[File]] input " should {
    "accept an Array[Array[String]] as the value for the input" in {
      runWdlAndAssertOutputs(
        sampleWdl = SampleWdl.ArrayOfArrays,
        eventFilter = EventFilter.info(pattern = "Workflow complete", occurrences = 1),
        expectedOutputs = Map(
          "wf.subtask.concatenated" -> WdlArray(WdlArrayType(WdlStringType), Seq(
            WdlString("foo\nbar\nbaz"),
            WdlString("third\nfourth")
          ))
        )
      )
    }
  }
}
