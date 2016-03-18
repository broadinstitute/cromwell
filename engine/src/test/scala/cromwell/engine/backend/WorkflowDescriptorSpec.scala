package cromwell.engine.backend

import java.nio.file.{Files, Paths}

import better.files._
import com.typesafe.config.ConfigFactory
import cromwell.CromwellTestkitSpec
import cromwell.core.WorkflowId
import cromwell.engine.backend.io._
import cromwell.engine.{ExecutionStatus, WorkflowSourceFiles, WorkflowSucceeded}
import cromwell.util.SampleWdl
import cromwell.webservice.WorkflowMetadataResponse
import org.joda.time.DateTime
import org.scalatest.concurrent.ScalaFutures
import org.scalatest.prop.TableDrivenPropertyChecks._
import org.scalatest.prop.Tables.Table
import org.scalatest.time.{Millis, Seconds, Span}
import org.scalatest.{FlatSpec, Matchers}
import spray.json.JsObject
import wdl4s.types.{WdlArrayType, WdlFileType}
import wdl4s.values.{WdlArray, WdlFile}

import scala.collection.JavaConverters._

class WorkflowDescriptorSpec extends FlatSpec with Matchers with ScalaFutures {

  val dummyWdl = "workflow w {}"
  val defaultConfig = Map("backend.backend" -> "local")
  val configWithCallCachingOn = defaultConfig + ("call-caching.enabled" -> "true")
  val configWithCallCachingOff = defaultConfig + ("call-caching.enabled" -> "false")
  val backend = new CromwellTestkitSpec.TestWorkflowManagerSystem()

  it should "honor configuration and workflow options for call-caching" in {
    val configs = Seq(defaultConfig, configWithCallCachingOn, configWithCallCachingOff)
    val options = Seq(None, Some(true), Some(false))

    val allCombinations = (for {
      config <- configs
      writeOption <- options
      readOption <- options
    } yield (config, writeOption, readOption)).toSet

    // writeCache is ON when config is ON and write-to-cache is None or true
    val writeCacheOnCombinations = (for {
      config <- configs if config == configWithCallCachingOn
      writeOption <- options if writeOption.isEmpty || writeOption.get
      readOption <- options
    } yield (config, writeOption, readOption)).toSet

    // readCache is ON when config is ON and read_from_cache is None or true
    val readCacheOnCombinations = (for {
      config <- configs if config == configWithCallCachingOn
      writeOption <- options
      readOption <- options if readOption.isEmpty || readOption.get
    } yield (config, writeOption, readOption)).toSet

    def makeOptions(writeOpt: Option[Boolean], readOpt: Option[Boolean]) = {
      val writeValue = writeOpt map { v => s""""write_to_cache": $v""" }
      val readValue = readOpt map { v => s""""read_from_cache": $v""" }
      val workflowOptions = Seq(writeValue, readValue).flatten.mkString(",\n")

      s"{$workflowOptions}"
    }

    def makeWorkflowDescriptor(config: Map[String, String], write: Option[Boolean], read: Option[Boolean]) = {
      val sources = WorkflowSourceFiles(dummyWdl, "{}", makeOptions(write, read))
      WorkflowDescriptor(WorkflowId.randomId(), sources, ConfigFactory.parseMap(config.asJava))
    }

    writeCacheOnCombinations foreach {
      case (config, writeToCache, readFromCache) => makeWorkflowDescriptor(config, writeToCache, readFromCache).writeToCache shouldBe true
    }

    readCacheOnCombinations foreach {
      case (config, writeToCache, readFromCache) => makeWorkflowDescriptor(config, writeToCache, readFromCache).readFromCache shouldBe true
    }

    (allCombinations -- writeCacheOnCombinations) foreach {
      case (config, writeToCache, readFromCache) => makeWorkflowDescriptor(config, writeToCache, readFromCache).writeToCache shouldBe false
    }

    (allCombinations -- readCacheOnCombinations) foreach {
      case (config, writeToCache, readFromCache) => makeWorkflowDescriptor(config, writeToCache, readFromCache).readFromCache shouldBe false
    }

  }

  def workflowFile(descriptor: WorkflowDescriptor, file: String): WdlFile = {
    val path = descriptor.wfContext.root.toPath(defaultFileSystems).resolve(file).toAbsolutePath
    path.createIfNotExists() // TODO: Delete the files and directories!
    WdlFile(path.toString)
  }

  it should "copy workflow outputs to their final (local) destination" in {
    import scala.concurrent.ExecutionContext.Implicits.global

    val tmpDir = Files.createTempDirectory("wf-outputs").toAbsolutePath
    val sources = WorkflowSourceFiles(SampleWdl.WorkflowOutputsWithFiles.wdlSource(), "{}",
      s"""{ "outputs_path": "$tmpDir" }""")

    val workflowOutputs = Map(
      "wfoutputs.A.out" -> "call-A/out",
      "wfoutputs.A.out2" -> "call-A/out2",
      "wfoutputs.B.outs" -> Seq("call-B/out", "call-B/out2"))

    val descriptor = WorkflowDescriptor(WorkflowId.randomId(), sources, ConfigFactory.load)

    val metadataOutputs = workflowOutputs mapValues {
      case file: String => workflowFile(descriptor, file)
      case files: Seq[_] => WdlArray(
        WdlArrayType(WdlFileType), files map {
          case file: String => workflowFile(descriptor, file)
        }
      )
    }

    val workflowMetadataResponse = WorkflowMetadataResponse(
      descriptor.id.toString,
      "wfoutputs",
      WorkflowSucceeded.toString,
      new DateTime(0),
      Option(new DateTime(0)),
      Option(new DateTime(0)),
      JsObject(),
      Option(metadataOutputs),
      Map.empty,
      None)

    implicit val patienceConfig = PatienceConfig(timeout = Span(10, Seconds), interval = Span(15, Millis))

    descriptor.copyWorkflowOutputs(workflowMetadataResponse).futureValue
    descriptor.maybeDeleteWorkflowLog()

    val expectedOutputs = Table(
      ("call", "file"),
      ("call-A", "out"),
      ("call-A", "out2"),
      ("call-B", "out"),
      ("call-B", "out2")
    )

    forAll(expectedOutputs) { (call, file) =>
      val path = tmpDir.resolve(Paths.get(descriptor.name, descriptor.id.toString, call, file))
      path.toFile should exist
    }
  }

  it should "copy call logs to their final (local) destination" in {
    import scala.concurrent.ExecutionContext.Implicits.global

    val tmpDir = Files.createTempDirectory("call-logs").toAbsolutePath
    val sources = WorkflowSourceFiles(SampleWdl.WorkflowOutputsWithFiles.wdlSource(), "{}",
      s"""{ "call_logs_dir": "$tmpDir" }""")
    val descriptor = WorkflowDescriptor(WorkflowId.randomId(), sources, ConfigFactory.load)

    val calls = Seq("call-A", "call-B").map({
      case call =>
        val metadata = CallMetadata(
          inputs=Map.empty,
          executionStatus=ExecutionStatus.Done.toString,
          backend=None,
          backendStatus=None,
          outputs=None,
          start=None,
          end=None,
          jobId=None,
          returnCode=None,
          shardIndex=0,
          stdout=Option(workflowFile(descriptor, s"$call/stdout")),
          stderr=Option(workflowFile(descriptor, s"$call/stderr")),
          backendLogs=if (call == "call-A") Option(Map("backendLog" -> workflowFile(descriptor, s"$call/backendout"))) else None,
          executionEvents=Seq.empty,
          attempt=1,
          runtimeAttributes=Map.empty,
          preemptible=None,
          cache=None,
          failures=None)
        call -> Seq(metadata)
    }).toMap

    val workflowMetadataResponse = WorkflowMetadataResponse(
      descriptor.id.toString,
      "wfoutputs",
      WorkflowSucceeded.toString,
      new DateTime(0),
      Option(new DateTime(0)),
      Option(new DateTime(0)),
      JsObject(),
      None,
      calls,
      None)

    descriptor.copyCallLogs(workflowMetadataResponse).futureValue
    descriptor.maybeDeleteWorkflowLog()

    val expectedOutputs = Table(
      ("call", "file", "checkExists"),
      ("call-A", "stdout", true),
      ("call-A", "stderr", true),
      ("call-A", "backendout", true),
      ("call-B", "stdout", true),
      ("call-B", "stderr", true),
      ("call-B", "backendout", false))

    forAll(expectedOutputs) { (call, file, checkExists) =>
      val path = tmpDir.resolve(Paths.get(descriptor.name, descriptor.id.toString, call, file))
      if (checkExists) {
        path.toFile should exist
      } else {
        path.toFile shouldNot exist
      }
    }
  }

  it should "copy workflow log to its final (local) destination" in {
    import scala.concurrent.ExecutionContext.Implicits.global

    val tmpDir = Files.createTempDirectory("workflow-log").toAbsolutePath
    val sources = WorkflowSourceFiles(SampleWdl.WorkflowOutputsWithFiles.wdlSource(), "{}",
      s"""{ "workflow_log_dir": "$tmpDir" }""")

    val descriptor = WorkflowDescriptor(WorkflowId.randomId(), sources, ConfigFactory.load)
    descriptor.copyWorkflowLog().futureValue
    descriptor.maybeDeleteWorkflowLog()
    val path = tmpDir.resolve(s"workflow.${descriptor.id}.log")
    path.toFile should exist
    descriptor.workflowLogPath shouldNot be(empty)
    descriptor.workflowLogPath.get.toFile shouldNot exist
  }

  it should "leave a workflow log in the temporary directory" in {
    import scala.concurrent.ExecutionContext.Implicits.global

    val tmpDir = Files.createTempDirectory("workflow-log").toAbsolutePath
    val sources = WorkflowSourceFiles(SampleWdl.WorkflowOutputsWithFiles.wdlSource(), "{}", "{}")

    val descriptor = WorkflowDescriptor(WorkflowId.randomId(), sources, ConfigFactory.parseString(
      s"""
         |workflow-options {
         |  workflow-log-dir: "$tmpDir"
         |  workflow-log-temporary: false
         |}""".stripMargin))
    descriptor.workflowLogPath shouldNot be(empty)
    descriptor.workflowLogPath.get.toFile shouldNot exist
    descriptor.workflowLogPath.foreach(_.createIfNotExists())
    descriptor.copyWorkflowLog().futureValue
    descriptor.maybeDeleteWorkflowLog()
    descriptor.workflowLogPath.get.toFile should exist
  }

  it should "not create a workflow log for an empty directory path" in {
    import scala.concurrent.ExecutionContext.Implicits.global

    val sources = WorkflowSourceFiles(SampleWdl.WorkflowOutputsWithFiles.wdlSource(), "{}", "{}")

    val descriptor = WorkflowDescriptor(WorkflowId.randomId(), sources, ConfigFactory.parseString(
      """workflow-options.workflow-log-dir: """""))
    descriptor.workflowLogPath should be(empty)
    descriptor.copyWorkflowLog().futureValue
    descriptor.maybeDeleteWorkflowLog()
    descriptor.workflowLogPath should be(empty)
  }

  it should "build the workflow root path" in {
    val sources = WorkflowSourceFiles(SampleWdl.WorkflowOutputsWithFiles.wdlSource(), "{}", "{}")
    val randomId: WorkflowId = WorkflowId.randomId()
    val descriptor = WorkflowDescriptor(randomId, sources, ConfigFactory.load)
    descriptor.workflowRootPathWithBaseRoot("/root") shouldBe Paths.get(s"/root/wfoutputs/$randomId")
  }

}
