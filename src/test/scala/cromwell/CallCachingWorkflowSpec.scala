package cromwell

import akka.testkit._
import cromwell.binding.types.{WdlIntegerType, WdlStringType, WdlArrayType}
import cromwell.binding.values.{WdlInteger, WdlFile, WdlArray, WdlString}
import cromwell.CromwellSpec.DockerTest
import cromwell.util.SampleWdl

import scala.language.postfixOps

class CallCachingWorkflowSpec extends CromwellTestkitSpec("CallCachingWorkflowSpec") {
  "A workflow which is run twice" should {
    "use cached calls on the second run" in {

      /** This workflow has two identical calls in it.
        *
        * Therefore the first invocation will get ONE call cache hit to itself
        * The second invocation will get TWO call cache hits.
        *
        * Altogether there are THREE call cache hits to the exact same call
        */
      runWdlAndAssertOutputs(
        sampleWdl = SampleWdl.FilePassingWorkflow,
        eventFilter = EventFilter.info(pattern = s"Call Caching: Cache hit. Using UUID\\(.{8}\\):a\\.*", occurrences = 1),
        expectedOutputs = Map(
          "file_passing.a.out" -> WdlFile("out"),
          "file_passing.a.out_interpolation" -> WdlFile("out"),
          "file_passing.a.contents" -> WdlString("foo bar baz"),
          "file_passing.b.out" -> WdlFile("out"),
          "file_passing.b.out_interpolation" -> WdlFile("out"),
          "file_passing.b.contents" -> WdlString("foo bar baz")
        )
      )
      runWdlAndAssertOutputs(
        sampleWdl = SampleWdl.FilePassingWorkflow,
        eventFilter = EventFilter.info(pattern = s"Call Caching: Cache hit. Using UUID\\(.{8}\\):a\\.*", occurrences = 2),
        expectedOutputs = Map(
          "file_passing.a.out" -> WdlFile("out"),
          "file_passing.a.out_interpolation" -> WdlFile("out"),
          "file_passing.a.contents" -> WdlString("foo bar baz"),
          "file_passing.b.out" -> WdlFile("out"),
          "file_passing.b.out_interpolation" -> WdlFile("out"),
          "file_passing.b.contents" -> WdlString("foo bar baz")
        )
      )
    }

    "use cached calls on the second run (docker)" taggedAs DockerTest in {
      runWdlAndAssertOutputs(
        sampleWdl = SampleWdl.FilePassingWorkflow,
        eventFilter = EventFilter.info(pattern = s"Call Caching: Cache hit. Using UUID\\(.{8}\\):a\\.*", occurrences = 1),
        runtime =
          """runtime {
            |  docker: "ubuntu:latest"
            |}
          """.stripMargin,
        expectedOutputs = Map(
          "file_passing.a.out" -> WdlFile("out"),
          "file_passing.a.out_interpolation" -> WdlFile("out"),
          "file_passing.a.contents" -> WdlString("foo bar baz"),
          "file_passing.b.out" -> WdlFile("out"),
          "file_passing.b.out_interpolation" -> WdlFile("out"),
          "file_passing.b.contents" -> WdlString("foo bar baz")
        )
      )
      runWdlAndAssertOutputs(
        sampleWdl = SampleWdl.FilePassingWorkflow,
        eventFilter = EventFilter.info(pattern = s"Call Caching: Cache hit. Using UUID\\(.{8}\\):a\\.*", occurrences = 2),
        runtime =
          """runtime {
            |  docker: "ubuntu:latest"
            |}
          """.stripMargin,
        expectedOutputs = Map(
          "file_passing.a.out" -> WdlFile("out"),
          "file_passing.a.out_interpolation" -> WdlFile("out"),
          "file_passing.a.contents" -> WdlString("foo bar baz"),
          "file_passing.b.out" -> WdlFile("out"),
          "file_passing.b.out_interpolation" -> WdlFile("out"),
          "file_passing.b.contents" -> WdlString("foo bar baz")
        )
      )
    }

    "use cached calls on the second run (scatter)" in {
      runWdlAndAssertOutputs(
        sampleWdl = new SampleWdl.ScatterWdl,
        eventFilter = EventFilter.info(pattern = s"starting calls: w.A", occurrences = 1),
        expectedOutputs = Map(
          "w.E.E_out" -> WdlArray(WdlArrayType(WdlIntegerType), Seq(9, 9, 9, 9, 9, 9).map(WdlInteger(_))),
          "w.C.C_out" -> WdlArray(WdlArrayType(WdlIntegerType), Seq(400, 500, 600, 800, 600, 500).map(WdlInteger(_))),
          "w.A.A_out" -> WdlArray(WdlArrayType(WdlStringType), Seq("jeff", "chris", "miguel", "thibault", "khalid", "scott").map(WdlString)),
          "w.D.D_out" -> WdlInteger(34),
          "w.B.B_out" -> WdlArray(WdlArrayType(WdlIntegerType), Seq(4, 5, 6, 8, 6, 5).map(WdlInteger(_)))
        )
      )
      runWdlAndAssertOutputs(
        sampleWdl = new SampleWdl.ScatterWdl,
        eventFilter = EventFilter.info(pattern = s"Call Caching: Cache hit. Using UUID\\(.{8}\\):A\\.*", occurrences = 1),
        expectedOutputs = Map(
          "w.E.E_out" -> WdlArray(WdlArrayType(WdlIntegerType), Seq(9, 9, 9, 9, 9, 9).map(WdlInteger(_))),
          "w.C.C_out" -> WdlArray(WdlArrayType(WdlIntegerType), Seq(400, 500, 600, 800, 600, 500).map(WdlInteger(_))),
          "w.A.A_out" -> WdlArray(WdlArrayType(WdlStringType), Seq("jeff", "chris", "miguel", "thibault", "khalid", "scott").map(WdlString)),
          "w.D.D_out" -> WdlInteger(34),
          "w.B.B_out" -> WdlArray(WdlArrayType(WdlIntegerType), Seq(4, 5, 6, 8, 6, 5).map(WdlInteger(_)))
        )
      )
    }
  }
}
