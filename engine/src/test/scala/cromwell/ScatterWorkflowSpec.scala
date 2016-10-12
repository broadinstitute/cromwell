package cromwell

import akka.testkit._
import cromwell.core.Tags.DockerTest
import wdl4s.types.{WdlArrayType, WdlFileType, WdlIntegerType, WdlStringType}
import wdl4s.values.{WdlArray, WdlFile, WdlInteger, WdlString}
import cromwell.util.SampleWdl

class ScatterWorkflowSpec extends CromwellTestKitSpec {
  "A workflow with a stand-alone scatter block in it" should {
    "run properly" in {
      runWdlAndAssertOutputs(
        sampleWdl = SampleWdl.SimpleScatterWdl,
        eventFilter = EventFilter.info(pattern = "Workflow complete", occurrences = 1),
        expectedOutputs = Map(
          "scatter0_outside_scatter_out" -> WdlInteger(8000),
          "scatter0_inside_scatter_out" -> WdlArray(WdlArrayType(WdlIntegerType), Seq(1, 2, 3, 4, 5).map(WdlInteger(_)))
        )
      )
    }
  }
  "A workflow with multiple calls in the scatter block" should {
    "run properly" in {
      runWdlAndAssertOutputs(
        sampleWdl = new SampleWdl.ScatterWdl,
        eventFilter = EventFilter.info(pattern = "Workflow complete", occurrences = 1),
        expectedOutputs = Map(
          "w_E_E_out" -> WdlArray(WdlArrayType(WdlIntegerType), Seq(9, 9, 9, 9, 9, 9).map(WdlInteger(_))),
          "w_C_C_out" -> WdlArray(WdlArrayType(WdlIntegerType), Seq(400, 500, 600, 800, 600, 500).map(WdlInteger(_))),
          "w_A_A_out" -> WdlArray(WdlArrayType(WdlStringType), Seq("jeff", "chris", "miguel", "thibault", "khalid", "ruchi").map(WdlString)),
          "w_D_D_out" -> WdlInteger(34),
          "w_B_B_out" -> WdlArray(WdlArrayType(WdlIntegerType), Seq(4, 5, 6, 8, 6, 5).map(WdlInteger(_)))
        )
      )
    }
  }
  "A workflow with sibling scatter blocks" should {
    "run properly" in {
      runWdlAndAssertOutputs(
        sampleWdl = SampleWdl.SiblingsScatterWdl,
        eventFilter = EventFilter.info(pattern = "Workflow complete", occurrences = 1),
        expectedOutputs = Map(
          "w_E_E_out" -> WdlArray(WdlArrayType(WdlIntegerType), Seq(9, 9, 9, 9, 9, 9).map(WdlInteger(_))),
          "w_F_B_out" -> WdlArray(WdlArrayType(WdlIntegerType), Seq(4, 5, 6, 8, 6, 5).map(WdlInteger(_))),
          "w_C_C_out" -> WdlArray(WdlArrayType(WdlIntegerType), Seq(400, 500, 600, 800, 600, 500).map(WdlInteger(_))),
          "w_A_A_out" -> WdlArray(WdlArrayType(WdlStringType), Seq("jeff", "chris", "miguel", "thibault", "khalid", "ruchi").map(WdlString)),
          "w_D_D_out" -> WdlInteger(34),
          "w_B_B_out" -> WdlArray(WdlArrayType(WdlIntegerType), Seq(4, 5, 6, 8, 6, 5).map(WdlInteger(_)))
        )
      )
    }
  }

  "A workflow with scatter blocks and File inputs/outputs" should {
    "run properly" in {
      runWdlAndAssertOutputs(
        sampleWdl = SampleWdl.PrepareScatterGatherWdl(),
        eventFilter = EventFilter.info(pattern = "Workflow complete", occurrences = 1),
        expectedOutputs = Map(
          "sc_test_do_gather_sum" -> WdlInteger(11),
          "sc_test_do_prepare_split_files" -> WdlArray(WdlArrayType(WdlFileType), Seq("temp_aa", "temp_ab", "temp_ac", "temp_ad").map(WdlFile(_))),
          "sc_test_do_scatter_count_file" -> WdlArray(WdlArrayType(WdlFileType), (1 to 4).map(_ => WdlFile("output.txt")))
        )
      )
    }

    "run properly in a Docker environment" taggedAs DockerTest in {
      runWdlAndAssertOutputs(
        sampleWdl = SampleWdl.PrepareScatterGatherWdl(),
        eventFilter = EventFilter.info(pattern = "Workflow complete", occurrences = 1),
        runtime = """
                  |runtime {
                  |  docker: "ubuntu:latest"
                  |}
                  """.stripMargin,
        expectedOutputs = Map(
          "sc_test_do_gather_sum" -> WdlInteger(11),
          "sc_test_do_prepare_split_files" -> WdlArray(WdlArrayType(WdlFileType), Seq("temp_aa", "temp_ab", "temp_ac", "temp_ad").map(WdlFile(_))),
          "sc_test_do_scatter_count_file" -> WdlArray(WdlArrayType(WdlFileType), (1 to 4).map(_ => WdlFile("output.txt")))
        )
      )
    }
  }
}
