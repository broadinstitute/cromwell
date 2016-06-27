package cromwell.backend.caching

import com.typesafe.config.ConfigFactory
import org.scalatest.{FlatSpec, Matchers}
import scala.collection.JavaConverters._

// TODO PBE Adapt to how new caching works, but the test logic should not need much change
class CachingConfigSpec extends FlatSpec with Matchers {

  val defaultConfig = Map("backend.backend" -> "local")
  val configWithCallCachingOn = defaultConfig + ("call-caching.enabled" -> "true")
  val configWithCallCachingOff = defaultConfig + ("call-caching.enabled" -> "false")

  "CachingConfigSpec" should "honor configuration and workflow options for caching" ignore {
    val configs = Seq(defaultConfig, configWithCallCachingOn, configWithCallCachingOff)
    val options = Seq(None, Some(true), Some(false))

    def makeOptions(writeOpt: Option[Boolean], readOpt: Option[Boolean]) = {
      val writeValue = writeOpt map { v => s""""write_to_cache": $v""" }
      val readValue = readOpt map { v => s""""read_from_cache": $v""" }
      val workflowOptions = Seq(writeValue, readValue).flatten.mkString(",\n")

      s"{$workflowOptions}"
    }

    val allCombinations = (for {
      config <- configs
      writeOption <- options
      readOption <- options
    } yield (ConfigFactory.parseMap(config.asJava), makeOptions(writeOption, readOption))).toSet

    // writeCache is ON when config is ON and write-to-cache is None or true
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

    writeCacheOnCombinations foreach {
      case (config, wfOptions) => //cacheConfig(config, wfOptions).writeToCache shouldBe true
    }

    readCacheOnCombinations foreach {
      case (config, wfOptions) => //cacheConfig(config, wfOptions).readFromCache shouldBe true
    }

    (allCombinations -- writeCacheOnCombinations) foreach {
      case (config, wfOptions) => //cacheConfig(config, wfOptions).writeToCache shouldBe false
    }

    (allCombinations -- readCacheOnCombinations) foreach {
      case (config, wfOptions) => //cacheConfig(config, wfOptions).readFromCache shouldBe false
    }
  }
}
