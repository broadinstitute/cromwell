package cromwell

import akka.testkit._
import cromwell.core.Tags.PostWomTest
import cromwell.util.SampleWdl
import wom.types._
import wom.values._


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
