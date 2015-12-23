package cromwell

import java.util.UUID

import akka.testkit._
import cromwell.binding.types.{WdlIntegerType, WdlStringType, WdlArrayType}
import cromwell.binding.values.{WdlInteger, WdlFile, WdlArray, WdlString}
import cromwell.CromwellSpec.DockerTest
import cromwell.util.SampleWdl

import scala.language.postfixOps

class CallCachingWorkflowSpec extends CromwellTestkitSpec {
  def cacheHitMessageForCall(name: String) = s"Call Caching: Cache hit. Using UUID\\(.{8}\\):$name\\.*"

  val expectedOutputs = Map(
    "file_passing.a.out" -> WdlFile("out"),
    "file_passing.a.out_interpolation" -> WdlFile("out"),
    "file_passing.a.contents" -> WdlString("foo bar baz"),
    "file_passing.b.out" -> WdlFile("out"),
    "file_passing.b.out_interpolation" -> WdlFile("out"),
    "file_passing.b.contents" -> WdlString("foo bar baz")
  )

  "A workflow which is run twice" should {
    "use cached calls on the second run" in {

      /** This workflow has two identical calls in it (run in serial)
        *
        * Therefore the first invocation will get ONE call cache hit to itself
        * The second invocation will get TWO call cache hits.
        *
        * Altogether there are THREE call cache hits to the exact same call
        */

      val salt = UUID.randomUUID().toString
      runWdlAndAssertOutputs(
        sampleWdl = SampleWdl.CallCachingWorkflow(salt),
        eventFilter = EventFilter.info(pattern = cacheHitMessageForCall("a"), occurrences = 1),
        expectedOutputs = expectedOutputs
      )
      runWdlAndAssertOutputs(
        sampleWdl = SampleWdl.CallCachingWorkflow(salt),
        eventFilter = EventFilter.info(pattern = cacheHitMessageForCall("a"), occurrences = 2),
        expectedOutputs = expectedOutputs
      )
    }

    "NOT use cached calls on the second run if read_from_cache workflow option is false" in {
      val salt = UUID.randomUUID().toString
      runWdlAndAssertOutputs(
        sampleWdl = SampleWdl.CallCachingWorkflow(salt),
        eventFilter = EventFilter.info(pattern = cacheHitMessageForCall("a"), occurrences = 0),
        workflowOptions = """{"read_from_cache": false}""",
        expectedOutputs = expectedOutputs
      )
      runWdlAndAssertOutputs(
        sampleWdl = SampleWdl.CallCachingWorkflow(salt),
        eventFilter = EventFilter.info(pattern = cacheHitMessageForCall("a"), occurrences = 0),
        workflowOptions = """{"read_from_cache": false}""",
        expectedOutputs = expectedOutputs
      )
    }

    "NOT use cached calls on the second run if write_to_cache workflow option is false" in {
      val salt = UUID.randomUUID().toString
      runWdlAndAssertOutputs(
        sampleWdl = SampleWdl.CallCachingWorkflow(salt),
        eventFilter = EventFilter.info(pattern = cacheHitMessageForCall("a"), occurrences = 0),
        workflowOptions = """{"write_to_cache": false}""",
        expectedOutputs = expectedOutputs
      )
      runWdlAndAssertOutputs(
        sampleWdl = SampleWdl.CallCachingWorkflow(salt),
        eventFilter = EventFilter.info(pattern = cacheHitMessageForCall("a"), occurrences = 0),
        workflowOptions = """{"write_to_cache": false}""",
        expectedOutputs = expectedOutputs
      )
    }

    "use cached calls on the second run (docker)" taggedAs DockerTest in {
      val salt = UUID.randomUUID().toString
      runWdlAndAssertOutputs(
        sampleWdl = SampleWdl.CallCachingWorkflow(salt),
        eventFilter = EventFilter.info(pattern = cacheHitMessageForCall("a"), occurrences = 1),
        runtime =
          """runtime {
            |  docker: "ubuntu:latest"
            |}
          """.stripMargin,
        expectedOutputs = expectedOutputs
      )
      runWdlAndAssertOutputs(
        sampleWdl = SampleWdl.CallCachingWorkflow(salt),
        eventFilter = EventFilter.info(pattern = cacheHitMessageForCall("a"), occurrences = 2),
        runtime =
          """runtime {
            |  docker: "ubuntu:latest"
            |}
          """.stripMargin,
        expectedOutputs = expectedOutputs
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
        // occurences = 20 because B, C, E are scattered 6 ways and A, D are not scattered
        eventFilter = EventFilter.info(pattern = cacheHitMessageForCall("[ABCDE]"), occurrences = 20),
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
