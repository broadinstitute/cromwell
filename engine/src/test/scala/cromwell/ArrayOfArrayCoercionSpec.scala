package cromwell

import akka.testkit._
import wdl4s.types.{WdlArrayType, WdlStringType}
import wdl4s.values.{WdlArray, WdlString}
import cromwell.util.SampleWdl

import scala.language.postfixOps

class ArrayOfArrayCoercionSpec extends CromwellTestkitSpec {
  "A workflow that has an Array[Array[File]] input " should {
    "accept an Array[Array[String]] as the value for the input" ignore {
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
