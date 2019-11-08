package cromwell.services.metadata.hybridcarbonite

import akka.actor.ActorSystem
import com.typesafe.config.ConfigFactory
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

  it should "respect valid freeze-scan config settings" in {
    val config = ConfigFactory.parseString(
      """{
        |   enabled = true
        |   bucket = "my_test_bucket"
        |   filesystems {
        |     gcs {
        |       auth = "application-default"
        |     }
        |   }
        |   freeze-scan {
        |     initial-interval = 1 second
        |     max-interval = 5 seconds
        |     multiplier = 1.2
        |   }
        |}
      """.stripMargin
    )

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
    val config = ConfigFactory.parseString(
      """{
        |   enabled = true
        |   bucket = "my_test_bucket"
        |   filesystems {
        |     gcs {
        |       auth = "application-default"
        |     }
        |   }
        |   freeze-scan {
        |     max-interval = 1 second
        |     initial-interval = 5 seconds
        |     multiplier = 1.2
        |   }
        |}
      """.stripMargin
    )

    val carboniteConfig = HybridCarboniteConfig.parseConfig(config)

    carboniteConfig match {
      case Left(e) =>
        e.head shouldBe "max-interval must be greater than or equal to initial-interval in Carboniter 'freeze-scan' stanza"
      case Right(_) => fail(s"Expected to fail but the config was parsed correctly!")
    }
  }

  it should "reject intervals that are not durations in freeze-scan config settings" in {
    val config = ConfigFactory.parseString(
      """{
        |   enabled = true
        |   bucket = "my_test_bucket"
        |   filesystems {
        |     gcs {
        |       auth = "application-default"
        |     }
        |   }
        |   freeze-scan {
        |     initial-interval = I
        |     max-interval = like
        |     multiplier = turtles
        |   }
        |}
      """.stripMargin
    )

    val carboniteConfig = HybridCarboniteConfig.parseConfig(config)

    carboniteConfig match {
      case Left(e) =>
        e.head shouldBe "Failed to parse Carboniter 'freeze-scan' stanza (reason 1 of 1): String: 10: Invalid value at 'initial-interval': No number in duration value 'I'"
      case Right(_) => fail(s"Expected to fail but the config was parsed correctly!")
    }
  }
}
