package cromwell.engine

import java.nio.file.{Files, Paths}

import com.typesafe.config.ConfigFactory
import cromwell.util.SampleWdl
import org.scalatest.prop.TableDrivenPropertyChecks._
import org.scalatest.prop.Tables.Table
import org.scalatest.{FlatSpec, Matchers}

import scala.collection.JavaConverters._
import scala.util.{Failure, Success}

class WorkflowDescriptorSpec extends FlatSpec with Matchers {

  val dummyWdl = "workflow w {}"
  val defaultConfig = Map("backend.backend" -> "local")
  val configWithCallCachingOn = defaultConfig + ("call-caching.enabled" -> "true")
  val configWithCallCachingOff = defaultConfig + ("call-caching.enabled" -> "false")

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

  it should "copy workflow outputs to their final (local) destination" in {
    import scala.concurrent.ExecutionContext.Implicits.global

    val tmpDir = Files.createTempDirectory("wf-outputs").toAbsolutePath
    val sources = WorkflowSourceFiles(SampleWdl.WorkflowOutputsWithFiles.wdlSource(), "{}", s"""{ "outputs_path": "$tmpDir" }""")

    val outputs = Table(
      ("call", "file"),
      ("call-A", "out"),
      ("call-A", "out2"),
      ("call-B", "out"),
      ("call-B", "out2")
    )

    val descriptor = WorkflowDescriptor(WorkflowId.randomId(), sources, ConfigFactory.load)

    descriptor.postProcessWorkflow onComplete {
      case Success(_) =>
        forAll(outputs) { (call, file) =>
          val path = tmpDir.resolve(Paths.get(descriptor.name, descriptor.id.toString, call, file))
          Files.exists(path) shouldBe true
        }
      case Failure(f) => fail(f)
    }
  }

}
