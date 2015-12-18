package cromwell

import akka.testkit._
import cromwell.binding.types.{WdlFileType, WdlIntegerType, WdlStringType, WdlArrayType}
import cromwell.binding.values.{WdlFile, WdlInteger, WdlArray, WdlString}
import cromwell.CromwellSpec.DockerTest
import cromwell.util.SampleWdl

import scala.language.postfixOps

class ReadTsvWorkflowSpec extends CromwellTestkitSpec("ReadTsvWorkflowSpec") {
  "A workflow with array/map indexes in expressions" should {
    "run locally" in {
      runWdlAndAssertOutputs(
        sampleWdl = SampleWdl.ReadTsvWdl,
        eventFilter = EventFilter.info(pattern = s"starting calls: test.output_file_table, test.output_matrix, test.output_table", occurrences = 1),
        expectedOutputs = Map(
          "test.output_matrix.matrix" -> WdlArray(
            WdlArrayType(WdlArrayType(WdlIntegerType)),
            Seq(
              WdlArray(WdlArrayType(WdlIntegerType), Seq(0, 1, 2).map(WdlInteger(_))),
              WdlArray(WdlArrayType(WdlIntegerType), Seq(3, 4, 5).map(WdlInteger(_))),
              WdlArray(WdlArrayType(WdlIntegerType), Seq(6, 7, 8).map(WdlInteger(_)))
            )
          ),
          "test.output_table.table" -> WdlArray(
            WdlArrayType(WdlArrayType(WdlStringType)),
            Seq(
              WdlArray(WdlArrayType(WdlStringType), Seq("col0", "col1", "col2").map(WdlString)),
              WdlArray(WdlArrayType(WdlStringType), Seq("a", "b", "c").map(WdlString)),
              WdlArray(WdlArrayType(WdlStringType), Seq("x", "y", "z").map(WdlString))
            )
          ),
          "test.output_file_table.table" -> WdlArray(
            WdlArrayType(WdlArrayType(WdlFileType)),
            Seq(
              WdlArray(WdlArrayType(WdlFileType), Seq("first", "second").map(WdlFile(_))),
              WdlArray(WdlArrayType(WdlFileType), Seq("third", "fourth").map(WdlFile(_)))
            )
          )
        )
      )
    }
  }
}
