package cromwell

import akka.testkit._
import cromwell.binding.types.{WdlFileType, WdlArrayType}
import cromwell.binding.values.{WdlInteger, WdlFile, WdlArray, WdlString}
import cromwell.util.SampleWdl

import scala.language.postfixOps

class ArrayWorkflowSpec extends CromwellTestkitSpec("ArrayWorkflowSpec") {
  "A task which contains a parameter " should {
    "accept an array for the value" in {
      runWdlAndAssertOutputs(
        sampleWdl = SampleWdl.ArrayIO,
        EventFilter.info(pattern = s"starting calls: wf.concat, wf.find", occurrences = 1),
        expectedOutputs = Map(
          "wf.count_lines.count" -> WdlInteger(3),
          "wf.count_lines_array.count" -> WdlInteger(3)
        )
      )
    }
  }
}
