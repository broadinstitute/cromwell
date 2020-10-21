package cromwell.filesystems.ftp

import cloud.nio.impl.ftp.{FtpAnonymousCredentials, FtpAuthenticatedCredentials}
import com.typesafe.config.ConfigFactory
import common.assertion.CromwellTimeoutSpec
import org.scalatest.flatspec.AnyFlatSpec
import org.scalatest.matchers.should.Matchers

class FtpInstanceConfigurationSpec extends AnyFlatSpec with CromwellTimeoutSpec with Matchers {

  behavior of "FtpConfigurationSpec"

  it should "parse anonymous credentials" in {
    FtpInstanceConfiguration(ConfigFactory.empty()).ftpCredentials shouldBe FtpAnonymousCredentials
  }

  it should "parse authenticated credentials" in {
    FtpInstanceConfiguration(ConfigFactory.parseString(
      """
        |auth {
        |  username = "me"
        |  password = "mot de passe"
        |}
      """.stripMargin)).ftpCredentials shouldBe FtpAuthenticatedCredentials("me", "mot de passe", None)

    FtpInstanceConfiguration(ConfigFactory.parseString(
      """
        |auth {
        |  username = "me"
        |  password = "mot de passe"
        |  account = "account"
        |}
      """.stripMargin)).ftpCredentials shouldBe FtpAuthenticatedCredentials("me", "mot de passe", Option("account"))
  }

}
