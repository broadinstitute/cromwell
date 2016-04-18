package cromwell.engine.backend

import java.nio.file.Paths

import akka.actor.ActorSystem
import better.files._
import com.typesafe.config.{Config, ConfigFactory}
import cromwell.CromwellTestkitSpec
import cromwell.core.WorkflowId
import cromwell.engine.WorkflowSourceFiles
import cromwell.engine.backend.io._
import cromwell.engine.workflow.MaterializeWorkflowDescriptorActor
import cromwell.engine.workflow.MaterializeWorkflowDescriptorActor.{MaterializeWorkflowDescriptorFailure, MaterializationResult, MaterializeWorkflowDescriptorSuccess}
import cromwell.util.{PromiseActor, SampleWdl}
import wdl4s.values.WdlFile

import scala.collection.JavaConverters._
import scala.concurrent.Await

trait WorkflowDescriptorBuilder {

  implicit val awaitTimeout = CromwellTestkitSpec.timeoutDuration
  implicit val actorSystem: ActorSystem

  def materializeWorkflowDescriptorFromSources(id: WorkflowId = WorkflowId.randomId(),
                                               workflowSources: WorkflowSourceFiles,
                                               conf: Config = ConfigFactory.load): WorkflowDescriptor = {
    import PromiseActor.EnhancedActorRef

    import scala.concurrent.ExecutionContext.Implicits.global

    val materializeWorkflowDescriptorActor = actorSystem.actorOf(MaterializeWorkflowDescriptorActor.props())

    val wfDesc = materializeWorkflowDescriptorActor.askNoTimeout(MaterializeWorkflowDescriptorActor.MaterializeWorkflow(id, workflowSources, conf)).
      mapTo[MaterializationResult]  map {
      case MaterializeWorkflowDescriptorSuccess(workflowDescriptor) => workflowDescriptor
      case MaterializeWorkflowDescriptorFailure(error) => throw error
    }
    val workflowDesc = Await.result(wfDesc, awaitTimeout)
    actorSystem.stop(materializeWorkflowDescriptorActor)
    workflowDesc
  }
}

class WorkflowDescriptorSpec extends CromwellTestkitSpec with WorkflowDescriptorBuilder {

  override implicit val actorSystem = system

  val dummyWdl = "workflow w {}"
  val defaultConfig = Map("backend.backend" -> "local")
  val configWithCallCachingOn = defaultConfig + ("call-caching.enabled" -> "true")
  val configWithCallCachingOff = defaultConfig + ("call-caching.enabled" -> "false")

  def workflowFile(descriptor: WorkflowDescriptor, file: String): WdlFile = {
    val path = descriptor.wfContext.root.toPath(defaultFileSystems).resolve(file).toAbsolutePath
    path.createIfNotExists() // TODO: Delete the files and directories!
    WdlFile(path.toString)
  }

  "WorkflowDescriptor" should {
    "honor configuration and workflow options for call-caching" in {
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
        materializeWorkflowDescriptorFromSources(workflowSources = sources, conf = ConfigFactory.parseMap(config.asJava))
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

    "build the workflow root path" in {
      val sources = WorkflowSourceFiles(SampleWdl.WorkflowOutputsWithFiles.wdlSource(), "{}", "{}")
      val randomId: WorkflowId = WorkflowId.randomId()
      val descriptor = materializeWorkflowDescriptorFromSources(id = randomId, workflowSources = sources, conf = ConfigFactory.load)
      descriptor.workflowRootPathWithBaseRoot("/root") shouldBe Paths.get(s"/root/wfoutputs/$randomId")
    }
  }
}
