package cromwell

import akka.testkit._
import cromwell.binding.types.{WdlArrayType, WdlStringType}
import cromwell.binding.values.{WdlArray, WdlString}
import cromwell.util.SampleWdl

import scala.language.postfixOps

class ArrayOfArrayCoercionSpec extends CromwellTestkitSpec("ArrayOfArrayCoercionSpec") {
  "A task which contains a parameter " should {
    "accept an array for the value" in {
      runWdlAndAssertOutputs(
        sampleWdl = SampleWdl.ArrayOfArrays,
        EventFilter.info(pattern = s"starting calls: wf.subtask, wf.subtask", occurrences = 1),
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
