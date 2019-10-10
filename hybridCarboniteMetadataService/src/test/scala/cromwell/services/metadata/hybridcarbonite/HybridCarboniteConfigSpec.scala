package cromwell.services.metadata.hybridcarbonite

import akka.actor.ActorSystem
import com.typesafe.config.ConfigFactory
import cromwell.core.TestKitSuite
import org.scalatest.{FlatSpecLike, Matchers}

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
      case Right(c) => {
        c.enabled shouldBe true
        c.bucket shouldBe "my_test_bucket"
        c.pathBuilders.head.name shouldBe "Google Cloud Storage"
      }
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
      case Right(c) => {
        c.enabled shouldBe false
        c.bucket shouldBe "my_test_bucket"
        c.pathBuilders.head.name shouldBe "Google Cloud Storage"
      }
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
}
