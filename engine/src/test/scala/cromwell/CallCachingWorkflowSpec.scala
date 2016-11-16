package cromwell

import java.util.UUID

import akka.testkit._
import com.typesafe.config.ConfigFactory
import cromwell.CallCachingWorkflowSpec._
import cromwell.core.Tags.{DockerTest, _}
import cromwell.util.SampleWdl
import wdl4s.types.{WdlArrayType, WdlIntegerType, WdlStringType}
import wdl4s.values.{WdlArray, WdlFile, WdlInteger, WdlString}


class CallCachingWorkflowSpec extends CromwellTestKitSpec {
  def cacheHitMessageForCall(name: String) = s"Call Caching: Cache hit. Using UUID\\(.{8}\\):$name\\.*"

  val expectedOutputs = Map(
    "file_passing.a.out" -> WdlFile("out"),
    "file_passing.a.out_interpolation" -> WdlFile("out"),
    "file_passing.a.contents" -> WdlString("foo bar baz"),
    "file_passing.a.stdoutContent" -> WdlArray(WdlArrayType(WdlStringType), Seq(WdlString("Something"))),
    "file_passing.b.out" -> WdlFile("out"),
    "file_passing.b.out_interpolation" -> WdlFile("out"),
    "file_passing.b.contents" -> WdlString("foo bar baz"),
    "file_passing.b.stdoutContent" -> WdlArray(WdlArrayType(WdlStringType), Seq(WdlString("Something")))
  )

  "A workflow which is run twice" should {
    "use cached calls on the second run" taggedAs PostMVP ignore {

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
        expectedOutputs = expectedOutputs,
        config = callCachingConfig
      )
      runWdlAndAssertOutputs(
        sampleWdl = SampleWdl.CallCachingWorkflow(salt),
        eventFilter = EventFilter.info(pattern = cacheHitMessageForCall("a"), occurrences = 2),
        expectedOutputs = expectedOutputs,
        config = callCachingConfig
      )
    }

    "NOT use cached calls on the second run if read_from_cache workflow option is false" taggedAs PostMVP ignore {
      val salt = UUID.randomUUID().toString
      runWdlAndAssertOutputs(
        sampleWdl = SampleWdl.CallCachingWorkflow(salt),
        eventFilter = EventFilter.info(pattern = cacheHitMessageForCall("a"), occurrences = 0),
        workflowOptions = """{"read_from_cache": false}""",
        expectedOutputs = expectedOutputs,
        config = callCachingConfig
      )
      runWdlAndAssertOutputs(
        sampleWdl = SampleWdl.CallCachingWorkflow(salt),
        eventFilter = EventFilter.info(pattern = cacheHitMessageForCall("a"), occurrences = 0),
        workflowOptions = """{"read_from_cache": false}""",
        expectedOutputs = expectedOutputs,
        config = callCachingConfig
      )
    }

    "NOT use cached calls on the second run if write_to_cache workflow option is false" taggedAs PostMVP ignore {
      val salt = UUID.randomUUID().toString
      runWdlAndAssertOutputs(
        sampleWdl = SampleWdl.CallCachingWorkflow(salt),
        eventFilter = EventFilter.info(pattern = cacheHitMessageForCall("a"), occurrences = 0),
        workflowOptions = """{"write_to_cache": false}""",
        expectedOutputs = expectedOutputs,
        config = callCachingConfig
      )
      runWdlAndAssertOutputs(
        sampleWdl = SampleWdl.CallCachingWorkflow(salt),
        eventFilter = EventFilter.info(pattern = cacheHitMessageForCall("a"), occurrences = 0),
        workflowOptions = """{"write_to_cache": false}""",
        expectedOutputs = expectedOutputs,
        config = callCachingConfig
      )
    }

    "use cached calls on the second run (docker)" taggedAs (DockerTest, PostMVP) ignore {
      val salt = UUID.randomUUID().toString
      runWdlAndAssertOutputs(
        sampleWdl = SampleWdl.CallCachingWorkflow(salt),
        eventFilter = EventFilter.info(pattern = cacheHitMessageForCall("a"), occurrences = 1),
        runtime =
          """runtime {
            |  docker: "ubuntu:latest"
            |}
          """.stripMargin,
        expectedOutputs = expectedOutputs,
        config = callCachingConfig
      )
      runWdlAndAssertOutputs(
        sampleWdl = SampleWdl.CallCachingWorkflow(salt),
        eventFilter = EventFilter.info(pattern = cacheHitMessageForCall("a"), occurrences = 2),
        runtime =
          """runtime {
            |  docker: "ubuntu:latest"
            |}
          """.stripMargin,
        expectedOutputs = expectedOutputs,
        config = callCachingConfig
      )
    }

    "use cached calls on the second run (scatter)" taggedAs PostMVP ignore {
      val outputs = Map(
        "w.E.E_out" -> WdlArray(WdlArrayType(WdlIntegerType), Seq(9, 9, 9, 9, 9, 9).map(WdlInteger(_))),
        "w.C.C_out" -> WdlArray(WdlArrayType(WdlIntegerType), Seq(400, 500, 600, 800, 600, 500).map(WdlInteger(_))),
        "w.A.A_out" -> WdlArray(WdlArrayType(WdlStringType), Seq("jeff", "chris", "miguel", "thibault", "khalid", "scott").map(WdlString)),
        "w.D.D_out" -> WdlInteger(34),
        "w.B.B_out" -> WdlArray(WdlArrayType(WdlIntegerType), Seq(4, 5, 6, 8, 6, 5).map(WdlInteger(_)))
      )

      runWdlAndAssertOutputs(
        sampleWdl = new SampleWdl.ScatterWdl,
        eventFilter = EventFilter.info(pattern = s"starting calls: w.A", occurrences = 1),
        expectedOutputs = outputs,
        config = callCachingConfig
      )
      runWdlAndAssertOutputs(
        sampleWdl = new SampleWdl.ScatterWdl,
        // occurences = 20 because B, C, E are scattered 6 ways and A, D are not scattered
        eventFilter = EventFilter.info(pattern = cacheHitMessageForCall("[ABCDE]"), occurrences = 20),
        expectedOutputs = outputs,
        config = callCachingConfig
      )
    }

    "show valid values for call caching in metadata" taggedAs PostMVP ignore {
      /*
        FIXME: This test had been constructing a custom WorkflowManagerActor. I don't believe this is still necessary
        but this test is being ignored so I'm not sure
       */
//      val workflowId = runWdlAndAssertOutputs(
//        sampleWdl = SampleWdl.CallCachingWorkflow(UUID.randomUUID().toString),
//        eventFilter = EventFilter.info(pattern = cacheHitMessageForCall("a"), occurrences = 1),
//        expectedOutputs = expectedOutputs,
//        config = CallCachingWorkflowSpec.callCachingConfig)

//      val status = messageAndWait[WorkflowManagerStatusSuccess](WorkflowStatus(workflowId)).state
//      status shouldEqual WorkflowSucceeded
//
//      val metadata = messageAndWait[WorkflowManagerWorkflowMetadataSuccess](WorkflowMetadata(workflowId)).response
//      metadata should not be null
//
//      val callA = metadata.calls.get("file_passing.a").get.head
//      val callB = metadata.calls.get("file_passing.b").get.head
//
//      callA.cache.get.allowResultReuse shouldEqual true
//      callA.cache.get.cacheHitWorkflow shouldEqual None
//      callA.cache.get.cacheHitCall shouldEqual None
//      callB.cache.get.allowResultReuse shouldEqual true
//      callB.cache.get.cacheHitWorkflow shouldEqual Some(workflowId.id.toString)
//      callB.cache.get.cacheHitCall shouldEqual Some("file_passing.a")
    }
  }
}

object CallCachingWorkflowSpec {
  val callCachingConfig = ConfigFactory.parseString(
    s"""
       |call-caching {
       |  enabled = true
       |  lookup-docker-hash = true
       |}
     """.stripMargin)
    //.withFallback(OldStyleWorkflowManagerActor.defaultConfig)
}
