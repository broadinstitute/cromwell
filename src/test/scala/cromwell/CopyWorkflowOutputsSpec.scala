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
        eventFilter = EventFilter.info(pattern = "transitioning from Running to Succeeded.", occurrences = 1),
        runtime = "",
        workflowOptions = s""" { "outputs_path": "$tmpDir" } """,
        expectedOutputs = Seq("A.out", "A.out2", "B.outs") map { o => ("wfoutputs." + o) ->  CromwellTestkitSpec.AnyValueIsFine } toMap,
        allowOtherOutputs = false
      )

      forAll(outputs) { (call, file) =>
        val path = tmpDir.resolve(Paths.get("wfoutputs", workflowId.id.toString, call, file))
        Files.exists(path) shouldBe true
      }
      val path = tmpDir.resolve(Paths.get("wfoutputs", workflowId.id.toString, "call-C", "out"))
      Files.exists(path) shouldBe false
    }
  }

}
