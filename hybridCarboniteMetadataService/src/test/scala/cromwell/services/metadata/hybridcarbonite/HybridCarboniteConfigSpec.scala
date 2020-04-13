package cromwell.services.metadata.hybridcarbonite

import akka.actor.ActorSystem
import com.typesafe.config.{Config, ConfigFactory}
import cromwell.core.TestKitSuite
import org.scalatest.{FlatSpecLike, Matchers}

import scala.concurrent.duration._
import scala.language.postfixOps

class HybridCarboniteConfigSpec extends TestKitSuite("HybridCarboniteConfigSpec") with FlatSpecLike with Matchers {

  behavior of "HybridCarboniteConfig"

  implicit val as: ActorSystem = system

  it should "parse a valid config correctly" in {
    val config = ConfigFactory.parseString(
      """{
        |   bucket = "my_test_bucket"
        |   filesystems {
        |     gcs {
        |       auth = "application-default"
        |     }
        |   }
        |   metadata-freezing {
        |     initial-interval: 5 seconds
        |   }
        |}
      """.stripMargin
    )

    val carboniteConfig = HybridCarboniteConfig.parseConfig(config)

    carboniteConfig match {
      case Left(e) => fail(s"Expected to parse correctly but got failure. Reason: $e")
      case Right(c) =>
        c.freezingConfig.enabled shouldBe true
        c.bucket shouldBe "my_test_bucket"
        c.pathBuilders.head.name shouldBe "Google Cloud Storage"
        //noinspection RedundantDefaultArgument
        val defaultFreezeScanConfig = MetadataFreezingConfig(
          initialInterval = 5 seconds,
          maxInterval = 5 minutes,
          multiplier = 1.1,
          minimumSummaryEntryId = None,
          debugLogging = true
        )
        c.freezingConfig shouldBe defaultFreezeScanConfig
    }
  }

  it should "parse correctly with freezing turned off if 'metadata-freezing.initial-interval' is not present in config" in {
    val config = ConfigFactory.parseString(
      """{
        |   bucket = "my_test_bucket"
        |   filesystems {
        |     gcs {
        |       auth = "application-default"
        |     }
        |   }
        |}
      """.stripMargin
    )

    val carboniteConfig = HybridCarboniteConfig.parseConfig(config)

    carboniteConfig match {
      case Left(e) => fail(s"Expected to parse correctly but got failure. Reason: $e")
      case Right(c) =>
        c.freezingConfig.enabled shouldBe false
        c.bucket shouldBe "my_test_bucket"
        c.pathBuilders.head.name shouldBe "Google Cloud Storage"
    }
  }

  it should "fail if 'bucket' is not mentioned" in {
    val config = ConfigFactory.parseString(
      """{
        |   filesystems {
        |     gcs {
        |       auth = "application-default"
        |     }
        |   }
        |   metadata-freezing {
        |     initial-interval: 5 seconds
        |   }
        |}
      """.stripMargin
    )

    val carboniteConfig = HybridCarboniteConfig.parseConfig(config)

    carboniteConfig match {
      case Left(e) => e.head shouldBe "Failed to parse Carboniter 'bucket' field from config (reason 1 of 1): No configuration setting found for key 'bucket'"
      case Right(_) => fail(s"Expected to fail but the config was parsed correctly!")
    }
  }

  it should "fail if gcs filesystem is not mentioned" in {
    val config = ConfigFactory.parseString(
      """{
        |   bucket = "my_test_bucket"
        |}
      """.stripMargin
    )

    val carboniteConfig = HybridCarboniteConfig.parseConfig(config)

    carboniteConfig match {
      case Left(e) => e.head shouldBe "Failed to parse Carboniter 'filesystems.gcs' field from config (reason 1 of 1): No configuration setting found for key 'filesystems'"
      case Right(_) => fail(s"Expected to fail but the config was parsed correctly!")
    }
  }

  it should "fail if gcs auth is not mentioned" in {
    val config = ConfigFactory.parseString(
      """{
        |   bucket = "my_test_bucket"
        |   filesystems.gcs {}
        |   metadata-freezing {
        |     initial-interval: 5 seconds
        |   }
        |}
      """.stripMargin
    )

    val carboniteConfig = HybridCarboniteConfig.parseConfig(config)

    carboniteConfig match {
      case Left(e) => e.head shouldBe "No configuration setting found for key 'auth'"
      case Right(_) => fail(s"Expected to fail but the config was parsed correctly!")
    }
  }

  private def buildFreezeScanConfig(initInterval: Any, maxInterval: Any, multiplier: Any): Config = {
    ConfigFactory.parseString(s"""{
      |   bucket = "my_test_bucket"
      |   filesystems {
      |     gcs {
      |       auth = "application-default"
      |     }
      |   }
      |   metadata-freezing {
      |     initial-interval = $initInterval
      |     max-interval = $maxInterval
      |     multiplier = $multiplier
      |   }
      |}
      """.stripMargin)
  }

  it should "respect valid metadata-freezing config settings" in {
    val config = buildFreezeScanConfig(initInterval = "1 second", maxInterval = "5 seconds", multiplier = "1.2")
    val carboniteConfig = HybridCarboniteConfig.parseConfig(config)

    carboniteConfig match {
      case Left(e) => fail(s"Expected to parse correctly but got failure. Reason: $e")
      case Right(c) =>
        val expectedConfig = MetadataFreezingConfig(
          initialInterval = 1 second,
          maxInterval = 5 seconds,
          multiplier = 1.2,
          minimumSummaryEntryId = None,
          debugLogging = true
        )
        c.freezingConfig shouldBe expectedConfig
    }
  }

  it should "reject max interval < initial interval in metadata-freezing config settings" in {
    val config = buildFreezeScanConfig(initInterval = "5 seconds", maxInterval = "1 seconds", multiplier = "1.2")
    val carboniteConfig = HybridCarboniteConfig.parseConfig(config)

    carboniteConfig match {
      case Left(e) =>
        e.head shouldBe "'max-interval' must be greater than or equal to a finite 'initial-interval' in Carboniter 'metadata-freezing' stanza"
      case Right(_) => fail(s"Expected to fail but the config was parsed correctly!")
    }
  }

  it should "reject invalid values in metadata-freezing config" in {
    val config = buildFreezeScanConfig(initInterval = "I", maxInterval = "like", multiplier = "turtles")
    val carboniteConfig = HybridCarboniteConfig.parseConfig(config)

    carboniteConfig match {
      case Left(e) =>
        e.size shouldBe 3
        e.toList.toSet shouldEqual Set(
          "Failed to parse Carboniter 'metadata-freezing' stanza (reason 1 of 3): format error I",
          "Failed to parse Carboniter 'metadata-freezing' stanza (reason 2 of 3): String: 10: Invalid value at 'max-interval': No number in duration value 'like'",
          "Failed to parse Carboniter 'metadata-freezing' stanza (reason 3 of 3): String: 11: multiplier has type STRING rather than NUMBER"
        )
      case Right(_) => fail(s"Expected to fail but the config was parsed correctly!")
    }
  }
}
