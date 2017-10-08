package cromwell

import akka.testkit._
import cromwell.core.Tags.PostWomTest
import wdl.types.{WdlArrayType, WdlStringType}
import wdl.values.{WdlArray, WdlString}
import cromwell.util.SampleWdl


class ArrayOfArrayCoercionSpec extends CromwellTestKitWordSpec {
  "A workflow that has an Array[Array[File]] input " should {
    // TODO WOM: Scatter workflow
    "accept an Array[Array[String]] as the value for the input" taggedAs PostWomTest ignore {
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
