package cromwell.services.metadata.hybridcarbonite

import akka.actor.ActorSystem
import cats.data.NonEmptyList
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
        |   enabled = true
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
        c.enabled shouldBe true
        c.bucket shouldBe "my_test_bucket"
        c.pathBuilders.head.name shouldBe "Google Cloud Storage"
        //noinspection RedundantDefaultArgument
        val defaultFreezeScanConfig = HybridCarboniteFreezeScanConfig(
          initialInterval = 5 seconds,
          maxInterval = 5 minutes,
          multiplier = 1.1
        )
        c.freezeScanConfig shouldBe defaultFreezeScanConfig
    }
  }

  it should "parse correctly if 'enabled' config is not present" in {
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
        c.enabled shouldBe false
        c.bucket shouldBe "my_test_bucket"
        c.pathBuilders.head.name shouldBe "Google Cloud Storage"
    }
  }

  it should "fail if 'bucket' is not mentioned" in {
    val config = ConfigFactory.parseString(
      """{
        |   enabled = true
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
      case Left(e) => e.head shouldBe "Failed to parse Carboniter 'bucket' field from config (reason 1 of 1): No configuration setting found for key 'bucket'"
      case Right(_) => fail(s"Expected to fail but the config was parsed correctly!")
    }
  }

  it should "fail if gcs filesystem is not mentioned" in {
    val config = ConfigFactory.parseString(
      """{
        |   enabled = true
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
        |   enabled = true
        |   bucket = "my_test_bucket"
        |   filesystems.gcs {}
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
      |   enabled = true
      |   bucket = "my_test_bucket"
      |   filesystems {
      |     gcs {
      |       auth = "application-default"
      |     }
      |   }
      |   freeze-scan {
      |     initial-interval = $initInterval
      |     max-interval = $maxInterval
      |     multiplier = $multiplier
      |   }
      |}
      """.stripMargin)
  }

  it should "respect valid freeze-scan config settings" in {
    val config = buildFreezeScanConfig(initInterval = "1 second", maxInterval = "5 seconds", multiplier = "1.2")
    val carboniteConfig = HybridCarboniteConfig.parseConfig(config)

    carboniteConfig match {
      case Left(e) => fail(s"Expected to parse correctly but got failure. Reason: $e")
      case Right(c) =>
        val expectedConfig = HybridCarboniteFreezeScanConfig(
          initialInterval = 1 second,
          maxInterval = 5 seconds,
          multiplier = 1.2
        )
        c.freezeScanConfig shouldBe expectedConfig
    }
  }

  it should "reject max interval < initial interval in freeze-scan config settings" in {
    val config = buildFreezeScanConfig(initInterval = "5 seconds", maxInterval = "1 seconds", multiplier = "1.2")
    val carboniteConfig = HybridCarboniteConfig.parseConfig(config)

    carboniteConfig match {
      case Left(e) =>
        e.head shouldBe "'max-interval' must be greater than or equal to 'initial-interval' in Carboniter 'freeze-scan' stanza"
      case Right(_) => fail(s"Expected to fail but the config was parsed correctly!")
    }
  }

  it should "reject invalid values in freeze-scan config" in {
    val config = buildFreezeScanConfig(initInterval = "I", maxInterval = "like", multiplier = "turtles")
    val carboniteConfig = HybridCarboniteConfig.parseConfig(config)

    carboniteConfig match {
      case Left(e) =>
        e.size shouldBe 3
        e shouldEqual NonEmptyList.of(
          "Failed to parse Carboniter 'freeze-scan' stanza (reason 1 of 3): String: 10: Invalid value at 'initial-interval': No number in duration value 'I'",
          "Failed to parse Carboniter 'freeze-scan' stanza (reason 2 of 3): String: 11: Invalid value at 'max-interval': No number in duration value 'like'",
          "Failed to parse Carboniter 'freeze-scan' stanza (reason 3 of 3): String: 12: multiplier has type STRING rather than NUMBER"
        )
      case Right(_) => fail(s"Expected to fail but the config was parsed correctly!")
    }
  }
}
