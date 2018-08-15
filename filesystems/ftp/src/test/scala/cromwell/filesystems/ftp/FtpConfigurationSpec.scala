package cromwell.filesystems.ftp

import cloud.nio.impl.ftp.{FtpAnonymousCredentials, FtpAuthenticatedCredentials}
import com.typesafe.config.ConfigFactory
import org.scalatest.{FlatSpec, Matchers}

class FtpConfigurationSpec extends FlatSpec with Matchers {

  behavior of "FtpConfigurationSpec"

  it should "parse anonymous credentials" in {
    FtpConfiguration(ConfigFactory.empty()).ftpCredentials shouldBe FtpAnonymousCredentials
  }

  it should "parse authenticated credentials" in {
    FtpConfiguration(ConfigFactory.parseString(
      """
        |auth {
        |  username = "me"
        |  password = "mot de passe"
        |}
      """.stripMargin)).ftpCredentials shouldBe FtpAuthenticatedCredentials("me", "mot de passe", None)

    FtpConfiguration(ConfigFactory.parseString(
      """
        |auth {
        |  username = "me"
        |  password = "mot de passe"
        |  account = "account"
        |}
      """.stripMargin)).ftpCredentials shouldBe FtpAuthenticatedCredentials("me", "mot de passe", Option("account"))
  }

}
