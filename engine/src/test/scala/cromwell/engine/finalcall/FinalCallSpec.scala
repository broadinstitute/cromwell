package cromwell.engine.finalcall

import java.nio.file.{Files, Path, Paths}

import akka.testkit.EventFilter
import better.files._
import com.typesafe.config.{Config, ConfigFactory}
import cromwell.CromwellTestkitSpec
import cromwell.engine.workflow.{MaterializeWorkflowDescriptorActor, WorkflowManagerActor}
import cromwell.util.SampleWdl
import org.scalatest.prop.TableDrivenPropertyChecks._
import org.scalatest.prop.Tables.Table

class FinalCallSpec extends CromwellTestkitSpec {
  "FinalCall" should {
    "copy workflow outputs to their final (local) destination" in {
      val tmpDir = Files.createTempDirectory("wf-outputs").toAbsolutePath
      val (workflowId, _) = runWdlAndReturnIdOutputs(
        sampleWdl = SampleWdl.WorkflowOutputsWithFiles,
        workflowOptions = s"""{ "outputs_path": "$tmpDir" }""",
        eventFilter = EventFilter.info(pattern = ": transitioning from Running to Succeeded.", occurrences = 1)
      )

      val expectedOutputs = Table(
        ("call", "file"),
        ("call-A", "out"),
        ("call-A", "out2"),
        ("call-B", "out"),
        ("call-B", "out2")
      )

      forAll(expectedOutputs) { (call, file) =>
        val path = tmpDir.resolve(Paths.get("wfoutputs", workflowId.toString, call, file))
        path.toFile should exist
      }
    }

    "copy call logs to their final (local) destination" in {
      val tmpDir = Files.createTempDirectory("call-logs").toAbsolutePath
      val (workflowId, _) = runWdlAndReturnIdOutputs(
        sampleWdl = SampleWdl.WorkflowOutputsWithFiles,
        workflowOptions = s"""{ "call_logs_dir": "$tmpDir" }""",
        eventFilter = EventFilter.info(pattern = ": transitioning from Running to Succeeded.", occurrences = 1)
      )

      val expectedOutputs = Table(
        ("call", "file"),
        ("call-A", "stdout"),
        ("call-A", "stderr"),
        ("call-B", "stdout"),
        ("call-B", "stderr"))

      forAll(expectedOutputs) { (call, file) =>
        val path = tmpDir.resolve(Paths.get("wfoutputs", workflowId.toString, call, file))
        path.toFile should exist
      }
    }

    "copy workflow log to its final (local) destination" in {
      val tmpDir = Files.createTempDirectory("workflow-log").toAbsolutePath
      val (workflowId, _) = runWdlAndReturnIdOutputs(
        sampleWdl = SampleWdl.WorkflowOutputsWithFiles,
        workflowOptions = s"""{ "workflow_log_dir": "$tmpDir" }""",
        eventFilter = EventFilter.info(pattern = ": transitioning from Running to Succeeded.", occurrences = 1)
      )

      val workflowLogName = s"workflow.$workflowId.log"
      val workflowLogPath = workflowLogPathForConfig(workflowLogName, ConfigFactory.load)
      val path = tmpDir.resolve(workflowLogName)
      path.toFile should exist
      workflowLogPath shouldNot be(empty)
      /*
      The workflow actor postStop deletes the file. We'd need to know the actor system was shut down, and we don't
      isolate the actor system.

      workflowLogPath.get.toFile shouldNot exist
       */
    }

    "leave a workflow log in the temporary directory" in {
      val tmpDir = Files.createTempDirectory("workflow-log").toAbsolutePath

      val config = ConfigFactory.parseString(
        s"""
           |workflow-options {
           |  workflow-log-dir: "$tmpDir"
           |  workflow-log-temporary: false
           |}""".stripMargin)

      val (workflowId, _) = runWdlAndReturnIdOutputs(
        sampleWdl = SampleWdl.WorkflowOutputsWithFiles,
        workflowOptions = "{}",
        eventFilter = EventFilter.info(pattern = ": transitioning from Running to Succeeded.", occurrences = 1),
        config = config.withFallback(WorkflowManagerActor.defaultConfig)
      )

      val workflowLogName = s"workflow.$workflowId.log"
      val workflowLogPath = workflowLogPathForConfig(workflowLogName, config)
      val path = tmpDir.resolve(workflowLogName)
      path.toFile should exist
      workflowLogPath.get.toFile should exist
    }

    "not create a workflow log for an empty directory path" in {
      val config = ConfigFactory.parseString("""workflow-options.workflow-log-dir: """"")

      val (workflowId, _) = runWdlAndReturnIdOutputs(
        sampleWdl = SampleWdl.WorkflowOutputsWithFiles,
        workflowOptions = "{}",
        eventFilter = EventFilter.info(pattern = ": transitioning from Running to Succeeded.", occurrences = 1),
        config = config.withFallback(WorkflowManagerActor.defaultConfig)
      )

      val workflowLogName = s"workflow.$workflowId.log"
      val workflowLogPath = workflowLogPathForConfig(workflowLogName, config)
      workflowLogPath should be(empty)
    }
  }

  private def workflowLogPathForConfig(workflowLogName: String, config: Config): Option[Path] = {
    val workflowLogOptions = MaterializeWorkflowDescriptorActor.workflowLogOptions(config)
    workflowLogOptions.map(_.dir.createDirectories() / workflowLogName).map(_.path)
  }
}
