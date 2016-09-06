package cromwell.engine.workflow.lifecycle

import com.typesafe.config.{Config, ConfigFactory}
import cromwell.core.WorkflowOptions
import cromwell.core.callcaching.CallCachingMode
import org.scalatest.{FlatSpec, Matchers}

import scala.collection.JavaConverters._
import scalaz.{Failure => ScalazFailure, Success => ScalazSuccess}
import scala.util.{Success, Try}

class CachingConfigSpec extends FlatSpec with Matchers {

  val defaultConfig = Map("backend.backend" -> "local")
  val configWithCallCachingOn = defaultConfig + ("call-caching.enabled" -> "true")
  val configWithCallCachingOff = defaultConfig + ("call-caching.enabled" -> "false")

  behavior of "Validating call caching config"

    val configs = Seq(defaultConfig, configWithCallCachingOn, configWithCallCachingOff)
    val options = Seq(None, Some(true), Some(false))

    def makeOptions(writeOpt: Option[Boolean], readOpt: Option[Boolean]) = {
      val writeValue = writeOpt map { v => s""""write_to_cache": $v""" }
      val readValue = readOpt map { v => s""""read_from_cache": $v""" }
      val workflowOptions = Seq(writeValue, readValue).flatten.mkString(",\n")

      WorkflowOptions.fromJsonString(s"{$workflowOptions}")
    }

    val allCombinations = (for {
      config <- configs
      writeOption <- options
      readOption <- options
    } yield (ConfigFactory.parseMap(config.asJava), makeOptions(writeOption, readOption))).toSet

    // writeCache is ON when config is ON and write_to_cache is None or true
    val writeCacheOnCombinations = (for {
      config <- configs if config == configWithCallCachingOn
      writeOption <- options if writeOption.isEmpty || writeOption.get
      readOption <- options
    } yield (ConfigFactory.parseMap(config.asJava), makeOptions(writeOption, readOption))).toSet

    // readCache is ON when config is ON and read_from_cache is None or true
    val readCacheOnCombinations = (for {
      config <- configs if config == configWithCallCachingOn
      writeOption <- options
      readOption <- options if readOption.isEmpty || readOption.get
    } yield (ConfigFactory.parseMap(config.asJava), makeOptions(writeOption, readOption))).toSet

    val writeCacheOffCombinations = allCombinations -- writeCacheOnCombinations
    val readCacheOffCombinations = allCombinations -- readCacheOnCombinations

    validateCallCachingMode("write cache on options", writeCacheOnCombinations) { mode => mode.writeToCache should be(true) }
    validateCallCachingMode("read cache on options", readCacheOnCombinations) { mode => mode.readFromCache should be(true) }
    validateCallCachingMode("write cache off options", writeCacheOffCombinations) { mode => mode.writeToCache should be(false) }
    validateCallCachingMode("read cache off options", readCacheOffCombinations) { mode => mode.readFromCache should be(false) }


  private def validateCallCachingMode(testName: String, combinations: Set[(Config, Try[WorkflowOptions])])(verificationFunction: CallCachingMode => Unit) = {
    it should s"correctly identify $testName" in {
      combinations foreach {
        case (config, Success(wfOptions)) =>
          MaterializeWorkflowDescriptorActor.validateCallCachingMode(wfOptions, config) match {
            case ScalazSuccess(activity) => verificationFunction(activity)
            case ScalazFailure(errors) =>
              val errorsList = errors.list.toList.mkString(", ")
              fail(s"Failure generating Call Config Mode: $errorsList")
          }
        case x => fail(s"Unexpected test tuple: $x")
      }
    }
  }
}
