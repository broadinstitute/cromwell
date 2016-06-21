package cromwell

import java.nio.file.{Files, Paths}

import akka.testkit.EventFilter
import cromwell.util.SampleWdl
import org.scalatest.prop.TableDrivenPropertyChecks._
import org.scalatest.prop.Tables.Table

import scala.language.postfixOps

class CopyWorkflowOutputsSpec extends CromwellTestkitSpec {

  "CopyWorkflowOutputsCall" should {
    "copy workflow outputs" in {
      val workflowOutputsPath = "copy-workflow-outputs"

      val tmpDir = Files.createTempDirectory(workflowOutputsPath).toAbsolutePath

      val outputs = Table(
        ("call", "file"),
        ("call-A", "out"),
        ("call-A", "out2"),
        ("call-B", "out"),
        ("call-B", "out2")
      )

      val workflowId = runWdlAndAssertOutputs(
        sampleWdl = SampleWdl.WorkflowOutputsWithFiles,
        eventFilter = EventFilter.info(
          pattern = "transition from FinalizingWorkflowState to WorkflowSucceededState", occurrences = 1),
        runtime = "",
        workflowOptions = s""" { "final_workflow_outputs_dir": "$tmpDir" } """,
        expectedOutputs = Seq("A.out", "A.out2", "B.outs") map { o => ("wfoutputs." + o) -> CromwellTestkitSpec.AnyValueIsFine } toMap,
        allowOtherOutputs = false
      )

      forAll(outputs) { (call, file) =>
        val path = tmpDir.resolve(Paths.get("wfoutputs", workflowId.id.toString, call, file))
        path.toFile should exist
      }
      val path = tmpDir.resolve(Paths.get("wfoutputs", workflowId.id.toString, "call-C", "out"))
      path.toFile shouldNot exist
    }

    "copy scattered workflow outputs" in {
      val workflowOutputsPath = "copy-workflow-outputs"

      val tmpDir = Files.createTempDirectory(workflowOutputsPath).toAbsolutePath

      val shards = 0 to 9
      val outputNames = List("B1", "B2")

      val outputTuples = for {
        shard <- shards
        output <- outputNames
      } yield ("call-A", s"shard-$shard/$output")

      val outputs = Table(("call", "file"), outputTuples: _*)

      val workflowId = runWdlAndAssertOutputs(
        sampleWdl = SampleWdl.WorkflowScatterOutputsWithFileArrays,
        eventFilter = EventFilter.info(
          pattern = "transition from FinalizingWorkflowState to WorkflowSucceededState", occurrences = 1),
        runtime = "",
        workflowOptions = s""" { "final_workflow_outputs_dir": "$tmpDir" } """,
        expectedOutputs = Map("wfoutputs.A.outs" -> CromwellTestkitSpec.AnyValueIsFine),
        allowOtherOutputs = false
      )

      forAll(outputs) { (call, file) =>
        val path = tmpDir.resolve(Paths.get("wfoutputs", workflowId.id.toString, call, file))
        Files.exists(path) shouldBe true
      }
    }
  }
}
