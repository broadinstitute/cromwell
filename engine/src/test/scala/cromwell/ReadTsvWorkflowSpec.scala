package cromwell

import akka.testkit.EventFilter
import cromwell.util.SampleWdl
import wdl4s.types.{WdlFileType, WdlStringType, WdlIntegerType, WdlArrayType}
import wdl4s.values.{WdlFile, WdlString, WdlInteger, WdlArray}

class ReadTsvWorkflowSpec extends CromwellTestkitSpec {
  "A workflow with array/map indexes in expressions" should {
    "run locally" in {
      runWdlAndAssertOutputs(
        sampleWdl = SampleWdl.ReadTsvWdl,
        eventFilter = EventFilter.info(
          pattern = "Starting calls: test.output_file_table:NA:1, test.output_matrix:NA:1, test.output_table:NA:1",
          occurrences = 1),
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
