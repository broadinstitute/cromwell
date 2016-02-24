package cromwell.engine

import com.typesafe.config.{Config, ConfigFactory}
import org.scalatest.prop.TableDrivenPropertyChecks._
import org.scalatest.prop.Tables.Table
import org.scalatest.{FlatSpec, Matchers}

import scala.collection.JavaConverters._

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
      case (config, writeToCache, readFromCache) => {
        val x = makeWorkflowDescriptor(config, writeToCache, readFromCache)
        x.writeToCache shouldBe false
      }
    }

    (allCombinations -- readCacheOnCombinations) foreach {
      case (config, writeToCache, readFromCache) => makeWorkflowDescriptor(config, writeToCache, readFromCache).readFromCache shouldBe false
    }

  }

}
