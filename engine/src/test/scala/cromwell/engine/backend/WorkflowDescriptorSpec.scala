package cromwell.engine.backend

import java.nio.file.{Files, Paths}
import java.time.OffsetDateTime

import akka.actor.{Actor, ActorRef, ActorSystem, Props}
import better.files._
import com.typesafe.config.{Config, ConfigFactory}
import cromwell.CromwellTestkitSpec
import cromwell.core.{ExecutionStatus, WorkflowId, WorkflowSucceeded}
import cromwell.engine.backend.io._

import cromwell.engine.workflow.lifecycle.MaterializeWorkflowDescriptorActor
import cromwell.engine.workflow.lifecycle.MaterializeWorkflowDescriptorActor.{MaterializeWorkflowDescriptorCommand, MaterializeWorkflowDescriptorFailureResponse, MaterializeWorkflowDescriptorSuccessResponse, WorkflowDescriptorMaterializationResult}
import cromwell.engine.{EngineWorkflowDescriptor, WorkflowSourceFiles}
import cromwell.util.{PromiseActor, SampleWdl}
import org.scalatest.prop.TableDrivenPropertyChecks._
import org.scalatest.prop.Tables.Table
import spray.json.JsObject
import wdl4s.types.{WdlArrayType, WdlFileType}
import wdl4s.values.{WdlArray, WdlFile}

import scala.collection.JavaConverters._
import scala.concurrent.Await

trait WorkflowDescriptorBuilder {

  implicit val awaitTimeout = CromwellTestkitSpec.timeoutDuration
  implicit val actorSystem: ActorSystem

//  def materializeWorkflowDescriptorFromSources(id: WorkflowId = WorkflowId.randomId(),
//                                               workflowSources: WorkflowSourceFiles,
//                                               conf: Config = ConfigFactory.load): OldStyleWorkflowDescriptor = {
//    import PromiseActor.EnhancedActorRef
//    import scala.concurrent.ExecutionContext.Implicits.global
//
//    val materializeWorkflowDescriptorActor = actorSystem.actorOf(OldStyleMaterializeWorkflowDescriptorActor.props())
//    val wfDesc = materializeWorkflowDescriptorActor.askNoTimeout(
//      OldStyleMaterializeWorkflowDescriptorActor.MaterializeWorkflow(id, workflowSources, conf)
//    ).mapTo[MaterializationResult] map {
//      case MaterializeWorkflowDescriptorSuccess(workflowDescriptor) => workflowDescriptor
//      case MaterializeWorkflowDescriptorFailure(error) => throw error
//    }
//    val workflowDesc = Await.result(wfDesc, awaitTimeout)
//    actorSystem.stop(materializeWorkflowDescriptorActor)
//    workflowDesc
//  }

  def createMaterializedEngineWorkflowDescriptor(id: WorkflowId, workflowSources: WorkflowSourceFiles): EngineWorkflowDescriptor = {
    import akka.pattern.ask
    implicit val timeout = akka.util.Timeout(awaitTimeout)
    implicit val ec = actorSystem.dispatcher

    val serviceRegistryIgnorer = actorSystem.actorOf(Props.empty)
    val actor = actorSystem.actorOf(MaterializeWorkflowDescriptorActor.props(serviceRegistryIgnorer), "MaterializeWorkflowDescriptorActor-" + id.id)
    val workflowDescriptorFuture = actor.ask(
      MaterializeWorkflowDescriptorCommand(id, workflowSources, ConfigFactory.load)
    ).mapTo[WorkflowDescriptorMaterializationResult]

    Await.result(workflowDescriptorFuture map {
      case MaterializeWorkflowDescriptorSuccessResponse(workflowDescriptor) => workflowDescriptor
      case MaterializeWorkflowDescriptorFailureResponse(reason) => throw reason
    }, awaitTimeout)
  }
}

class WorkflowDescriptorSpec extends CromwellTestkitSpec with WorkflowDescriptorBuilder {

  override implicit val actorSystem = system

  val dummyWdl = "workflow w {}"
  val defaultConfig = Map("backend.backend" -> "local")
  val configWithCallCachingOn = defaultConfig + ("call-caching.enabled" -> "true")
  val configWithCallCachingOff = defaultConfig + ("call-caching.enabled" -> "false")

//  def workflowFile(descriptor: OldStyleWorkflowDescriptor, file: String): WdlFile = {
//    val path = descriptor.wfContext.root.toPath(defaultFileSystems).resolve(file).toAbsolutePath
//    path.createIfNotExists() // TODO: Delete the files and directories!
//    WdlFile(path.toString)
//  }

  "WorkflowDescriptor" should {
    "honor configuration and workflow options for call-caching" ignore {
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

//      def makeWorkflowDescriptor(config: Map[String, String], write: Option[Boolean], read: Option[Boolean]) = {
//        val sources = WorkflowSourceFiles(dummyWdl, "{}", makeOptions(write, read))
//        materializeWorkflowDescriptorFromSources(workflowSources = sources, conf = ConfigFactory.parseMap(config.asJava))
//      }

//      writeCacheOnCombinations foreach {
//        case (config, writeToCache, readFromCache) => makeWorkflowDescriptor(config, writeToCache, readFromCache).writeToCache shouldBe true
//      }
//
//      readCacheOnCombinations foreach {
//        case (config, writeToCache, readFromCache) => makeWorkflowDescriptor(config, writeToCache, readFromCache).readFromCache shouldBe true
//      }
//
//      (allCombinations -- writeCacheOnCombinations) foreach {
//        case (config, writeToCache, readFromCache) => makeWorkflowDescriptor(config, writeToCache, readFromCache).writeToCache shouldBe false
//      }
//
//      (allCombinations -- readCacheOnCombinations) foreach {
//        case (config, writeToCache, readFromCache) => makeWorkflowDescriptor(config, writeToCache, readFromCache).readFromCache shouldBe false
//      }

    }

    "copy workflow outputs to their final (local) destination" ignore {
//
//      val tmpDir = Files.createTempDirectory("wf-outputs").toAbsolutePath
//      val sources = WorkflowSourceFiles(SampleWdl.WorkflowOutputsWithFiles.wdlSource(), "{}",
//        s"""{ "outputs_path": "$tmpDir" }""")
//
//      val workflowOutputs = Map(
//        "wfoutputs.A.out" -> "call-A/out",
//        "wfoutputs.A.out2" -> "call-A/out2",
//        "wfoutputs.B.outs" -> Seq("call-B/out", "call-B/out2"))
//
//      val descriptor = materializeWorkflowDescriptorFromSources(workflowSources = sources, conf = ConfigFactory.load)
//
//      val metadataOutputs = workflowOutputs mapValues {
//        case file: String => workflowFile(descriptor, file)
//        case files: Seq[_] => WdlArray(
//          WdlArrayType(WdlFileType), files map {
//            case file: String => workflowFile(descriptor, file)
//          }
//        )
//      }
//
//      val workflowMetadataResponse = WorkflowMetadataResponse(
//        descriptor.id.toString,
//        "wfoutputs",
//        WorkflowSucceeded.toString,
//        OffsetDateTime.now,
//        Option(OffsetDateTime.now),
//        Option(OffsetDateTime.now),
//        JsObject(),
//        Option(metadataOutputs),
//        Map.empty)
//
//      descriptor.copyWorkflowOutputs(workflowMetadataResponse).futureValue
//      descriptor.maybeDeleteWorkflowLog()
//
//      val expectedOutputs = Table(
//        ("call", "file"),
//        ("call-A", "out"),
//        ("call-A", "out2"),
//        ("call-B", "out"),
//        ("call-B", "out2")
//      )
//
//      forAll(expectedOutputs) { (call, file) =>
//        val path = tmpDir.resolve(Paths.get(descriptor.name, descriptor.id.toString, call, file))
//        path.toFile should exist
//      }
    }

    "copy call logs to their final (local) destination" ignore {
//
//      val tmpDir = Files.createTempDirectory("call-logs").toAbsolutePath
//      val sources = WorkflowSourceFiles(SampleWdl.WorkflowOutputsWithFiles.wdlSource(), "{}",
//        s"""{ "call_logs_dir": "$tmpDir" }""")
//      val descriptor = materializeWorkflowDescriptorFromSources(workflowSources = sources, conf = ConfigFactory.load)
//
//    val calls = Seq("call-A", "call-B").map({
//      case call =>
//        val metadata = OldStyleCallMetadata(
//          inputs=Map.empty,
//          executionStatus=ExecutionStatus.Done.toString,
//          backend=None,
//          backendStatus=None,
//          outputs=None,
//          start=None,
//          end=None,
//          jobId=None,
//          returnCode=None,
//          shardIndex=0,
//          stdout=Option(workflowFile(descriptor, s"$call/stdout")),
//          stderr=Option(workflowFile(descriptor, s"$call/stderr")),
//          backendLogs=if (call == "call-A") Option(Map("backendLog" -> workflowFile(descriptor, s"$call/backendout"))) else None,
//          attempt=1,
//          runtimeAttributes=Map.empty,
//          preemptible=None,
//          cache=None,
//          failures=None)
//        call -> Seq(metadata)
//    }).toMap
//
//      val workflowMetadataResponse = WorkflowMetadataResponse(
//        descriptor.id.toString,
//        "wfoutputs",
//        WorkflowSucceeded.toString,
//        OffsetDateTime.now,
//        Option(OffsetDateTime.now),
//        Option(OffsetDateTime.now),
//        JsObject(),
//        None,
//        calls)
//
//      descriptor.copyCallLogs(workflowMetadataResponse).futureValue
//      descriptor.maybeDeleteWorkflowLog()
//
//      val expectedOutputs = Table(
//        ("call", "file", "checkExists"),
//        ("call-A", "stdout", true),
//        ("call-A", "stderr", true),
//        ("call-A", "backendout", true),
//        ("call-B", "stdout", true),
//        ("call-B", "stderr", true),
//        ("call-B", "backendout", false))
//
//      forAll(expectedOutputs) { (call, file, checkExists) =>
//        val path = tmpDir.resolve(Paths.get(descriptor.name, descriptor.id.toString, call, file))
//        if (checkExists) {
//          path.toFile should exist
//        } else {
//          path.toFile shouldNot exist
//        }
//      }
    }

    "copy workflow log to its final (local) destination" ignore {
//
//      val tmpDir = Files.createTempDirectory("workflow-log").toAbsolutePath
//      val sources = WorkflowSourceFiles(SampleWdl.WorkflowOutputsWithFiles.wdlSource(), "{}",
//        s"""{ "workflow_log_dir": "$tmpDir" }""")
//
//      val descriptor = materializeWorkflowDescriptorFromSources(workflowSources = sources, conf = ConfigFactory.load)
//      descriptor.copyWorkflowLog().futureValue
//      descriptor.maybeDeleteWorkflowLog()
//      val path = tmpDir.resolve(s"workflow.${descriptor.id}.log")
//      path.toFile should exist
//      descriptor.workflowLogPath shouldNot be(empty)
//      descriptor.workflowLogPath.get.toFile shouldNot exist
    }

    "leave a workflow log in the temporary directory" ignore {
//
//      val tmpDir = Files.createTempDirectory("workflow-log").toAbsolutePath
//      val sources = WorkflowSourceFiles(SampleWdl.WorkflowOutputsWithFiles.wdlSource(), "{}", "{}")
//
//      val descriptor = materializeWorkflowDescriptorFromSources(workflowSources = sources,
//        conf = ConfigFactory.parseString(
//          s"""
//             |workflow-options {
//             |  workflow-log-dir: "$tmpDir"
//             |  workflow-log-temporary: false
//             |}""".stripMargin))
//      descriptor.workflowLogPath shouldNot be(empty)
//      descriptor.workflowLogPath.get.toFile shouldNot exist
//      descriptor.workflowLogPath.foreach(_.createIfNotExists())
//      descriptor.copyWorkflowLog().futureValue
//      descriptor.maybeDeleteWorkflowLog()
//      descriptor.workflowLogPath.get.toFile should exist
    }

    "not create a workflow log for an empty directory path" ignore {
//
//      val sources = WorkflowSourceFiles(SampleWdl.WorkflowOutputsWithFiles.wdlSource(), "{}", "{}")
//
//      val descriptor = materializeWorkflowDescriptorFromSources(workflowSources = sources,
//        conf = ConfigFactory.parseString("""workflow-options.workflow-log-dir: """""))
//      descriptor.workflowLogPath should be(empty)
//      descriptor.copyWorkflowLog().futureValue
//      descriptor.maybeDeleteWorkflowLog()
//      descriptor.workflowLogPath should be(empty)
    }

    "build the workflow root path" ignore {
//      val sources = WorkflowSourceFiles(SampleWdl.WorkflowOutputsWithFiles.wdlSource(), "{}", "{}")
//      val randomId: WorkflowId = WorkflowId.randomId()
//      val descriptor = materializeWorkflowDescriptorFromSources(id = randomId, workflowSources = sources, conf = ConfigFactory.load)
//      descriptor.workflowRootPathWithBaseRoot("/root") shouldBe Paths.get(s"/root/wfoutputs/$randomId")
    }
  }
}
