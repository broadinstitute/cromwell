package cromwell

import akka.testkit._
import cromwell.core.Tags.{DockerTest, PostWomTest}
import cromwell.util.SampleWdl
import wom.types._
import wom.values._

// TODO WOM: need scatter support
class ScatterWorkflowSpec extends CromwellTestKitWordSpec {
  "A workflow with a stand-alone scatter block in it" should {
    "run properly" taggedAs PostWomTest ignore {
      runWdlAndAssertOutputs(
        sampleWdl = SampleWdl.SimpleScatterWdl,
        eventFilter = EventFilter.info(pattern = "Workflow complete", occurrences = 1),
        expectedOutputs = Map(
          "scatter0.outside_scatter.out" -> WomInteger(8000),
          "scatter0.inside_scatter.out" -> WomArray(WomArrayType(WomIntegerType), Seq(1, 2, 3, 4, 5).map(WomInteger(_)))
        )
      )
    }
  }
  "A workflow with multiple calls in the scatter block" should {
    "run properly" taggedAs PostWomTest ignore {
      runWdlAndAssertOutputs(
        sampleWdl = new SampleWdl.ScatterWdl,
        eventFilter = EventFilter.info(pattern = "Workflow complete", occurrences = 1),
        expectedOutputs = Map(
          "w.E.E_out" -> WomArray(WomArrayType(WomIntegerType), Seq(9, 9, 9, 9, 9, 9).map(WomInteger(_))),
          "w.C.C_out" -> WomArray(WomArrayType(WomIntegerType), Seq(400, 500, 600, 800, 600, 500).map(WomInteger(_))),
          "w.A.A_out" -> WomArray(WomArrayType(WomStringType), Seq("jeff", "chris", "miguel", "thibault", "khalid", "ruchi").map(WomString)),
          "w.D.D_out" -> WomInteger(34),
          "w.B.B_out" -> WomArray(WomArrayType(WomIntegerType), Seq(4, 5, 6, 8, 6, 5).map(WomInteger(_)))
        )
      )
    }
  }
  "A workflow with sibling scatter blocks" should {
    "run properly" taggedAs PostWomTest ignore {
      runWdlAndAssertOutputs(
        sampleWdl = SampleWdl.SiblingsScatterWdl,
        eventFilter = EventFilter.info(pattern = "Workflow complete", occurrences = 1),
        expectedOutputs = Map(
          "w.E.E_out" -> WomArray(WomArrayType(WomIntegerType), Seq(9, 9, 9, 9, 9, 9).map(WomInteger(_))),
          "w.F.B_out" -> WomArray(WomArrayType(WomIntegerType), Seq(4, 5, 6, 8, 6, 5).map(WomInteger(_))),
          "w.C.C_out" -> WomArray(WomArrayType(WomIntegerType), Seq(400, 500, 600, 800, 600, 500).map(WomInteger(_))),
          "w.A.A_out" -> WomArray(WomArrayType(WomStringType), Seq("jeff", "chris", "miguel", "thibault", "khalid", "ruchi").map(WomString)),
          "w.D.D_out" -> WomInteger(34),
          "w.B.B_out" -> WomArray(WomArrayType(WomIntegerType), Seq(4, 5, 6, 8, 6, 5).map(WomInteger(_)))
        )
      )
    }
  }

  "A workflow with scatter blocks and File inputs/outputs" should {
    "run properly" taggedAs PostWomTest ignore {
      runWdlAndAssertOutputs(
        sampleWdl = SampleWdl.PrepareScatterGatherWdl(),
        eventFilter = EventFilter.info(pattern = "Workflow complete", occurrences = 1),
        expectedOutputs = Map(
          "sc_test.do_gather.sum" -> WomInteger(11),
          "sc_test.do_prepare.split_files" -> WomArray(WomArrayType(WomFileType), Seq("temp_aa", "temp_ab", "temp_ac", "temp_ad").map(WomFile(_))),
          "sc_test.do_scatter.count_file" -> WomArray(WomArrayType(WomFileType), (1 to 4).map(_ => WomFile("output.txt")))
        )
      )
    }

    "run properly in a Docker environment" taggedAs (DockerTest, PostWomTest) ignore {
      runWdlAndAssertOutputs(
        sampleWdl = SampleWdl.PrepareScatterGatherWdl(),
        eventFilter = EventFilter.info(pattern = "Workflow complete", occurrences = 1),
        runtime = """
                  |runtime {
                  |  docker: "ubuntu:latest"
                  |}
                  """.stripMargin,
        expectedOutputs = Map(
          "sc_test.do_gather.sum" -> WomInteger(11),
          "sc_test.do_prepare.split_files" -> WomArray(WomArrayType(WomFileType), Seq("temp_aa", "temp_ab", "temp_ac", "temp_ad").map(WomFile(_))),
          "sc_test.do_scatter.count_file" -> WomArray(WomArrayType(WomFileType), (1 to 4).map(_ => WomFile("output.txt")))
        )
      )
    }
  }
}
