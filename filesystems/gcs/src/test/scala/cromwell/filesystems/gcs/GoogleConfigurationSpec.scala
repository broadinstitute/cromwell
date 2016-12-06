package cromwell.filesystems.gcs

import better.files.File
import com.typesafe.config.{ConfigException, ConfigFactory}
import cromwell.filesystems.gcs.GoogleConfiguration.GoogleConfigurationException
import cromwell.filesystems.gcs.auth.{ApplicationDefaultMode, RefreshTokenMode, ServiceAccountMode, UserMode}
import org.scalatest.{FlatSpec, Matchers}


class GoogleConfigurationSpec extends FlatSpec with Matchers {

  behavior of "GoogleConfiguration"

  it should "parse all manner of well-formed auths" in {
    val mockFile = File.newTemporaryFile()

    val righteousGoogleConfig =
      s"""
        |google {
        |  application-name = "cromwell"
        |
        |  auths = [
        |    {
        |      name = "name-default"
        |      scheme = "application_default"
        |    },
        |    {
        |      name = "name-refresh"
        |      scheme = "refresh_token"
        |      client-id = "secret_id"
        |      client-secret = "secret_secret"
        |    },
        |    {
        |      name = "name-user"
        |      scheme = "user_account"
        |      user = "me"
        |      secrets-file = "${mockFile.pathAsString}"
        |      data-store-dir = "/where/the/data/at"
        |    },
        |    {
        |      name = "name-service"
        |      scheme = "service_account"
        |      service-account-id = "my-google-account"
        |      pem-file = "${mockFile.pathAsString}"
        |    }
        |  ]
        |}
        |
      """.stripMargin

    val gconf = GoogleConfiguration(ConfigFactory.parseString(righteousGoogleConfig))

    gconf.applicationName shouldBe "cromwell"
    gconf.authsByName should have size 4

    val auths = gconf.authsByName.values

    val appDefault = (auths collectFirst { case a: ApplicationDefaultMode => a }).get
    appDefault.name shouldBe "name-default"

    val refreshToken = (auths collectFirst { case a: RefreshTokenMode => a }).get
    refreshToken.name shouldBe "name-refresh"
    refreshToken.clientSecret shouldBe "secret_secret"
    refreshToken.clientId shouldBe "secret_id"

    val user = (auths collectFirst { case a: UserMode => a }).get
    user.name shouldBe "name-user"
    user.secretsPath shouldBe mockFile.pathAsString
    user.datastoreDir shouldBe "/where/the/data/at"

    val service = (auths collectFirst { case a: ServiceAccountMode => a }).get
    service.name shouldBe "name-service"
    service.accountId shouldBe "my-google-account"
    service.pemPath shouldBe mockFile.pathAsString

    mockFile.delete(true)
  }


  it should "not parse a configuration stanza without applicationName" in {
    val applessGoogleConfig =
      """
        |google {
        |  auths = [
        |    {
        |      name = "name-default"
        |      scheme = "application_default"
        |    }
        |  ]
        |}
      """.stripMargin

    a[GoogleConfigurationException] shouldBe thrownBy {
      GoogleConfiguration(ConfigFactory.parseString(applessGoogleConfig))
    }
  }

  it should "not parse a configuration stanza with wrong cromwell auth" in {
    val unsupported =
      """
        |google {
        |  application-name = "cromwell"
        |
        |  auths = [
        |    {
        |      name = "unsupported-auth"
        |      scheme = "not supported"
        |    }
        |  ]
        |}
      """.stripMargin

    a[GoogleConfigurationException] shouldBe thrownBy {
      GoogleConfiguration(ConfigFactory.parseString(unsupported))
    }

    val schemeless =
      """
        |google {
        |  application-name = "cromwell"
        |
        |  auths = [
        |    {
        |      name = "scheme-unspecified"
        |    }
        |  ]
        |}
      """.stripMargin

    a[ConfigException.Missing] shouldBe thrownBy {
      GoogleConfiguration(ConfigFactory.parseString(schemeless))
    }

    val nameless =
      """
        |google {
        |  application-name = "cromwell"
        |
        |  auths = [
        |    {
        |      scheme = "application-default"
        |    }
        |  ]
        |}
      """.stripMargin

    a[ConfigException.Missing] shouldBe thrownBy {
      GoogleConfiguration(ConfigFactory.parseString(nameless))
    }

    // The various GoogleAuthModes actually don't complain about spurious keys in their
    // configurations as long as all the keys they do care about are present.  That's not
    // necessarily ideal behavior.
    val badKeyInRefreshTokenMode =
      """
        |google {
        |  application-name = "cromwell"
        |
        |  auths = [
        |    {
        |      name = "name-refresh"
        |      scheme = "refresh_token"
        |      client-id-botched-key = "secret_id"
        |      client-secret = "secret_secret"
        |    }
        |  ]
        |}
      """.stripMargin

    a[GoogleConfigurationException] shouldBe thrownBy {
      GoogleConfiguration(ConfigFactory.parseString(badKeyInRefreshTokenMode))
    }

    val badKeyInUserMode =
      """
        |google {
        |  application-name = "cromwell"
        |
        |  auths = [
        |    {
        |      name = "name-user"
        |      scheme = "user_account"
        |      user-botched-key = "me"
        |      secrets-file = "/very/secret/file.txt"
        |      data-store-dir = "/where/the/data/at"
        |    }
        |  ]
        |}
      """.stripMargin

    a[GoogleConfigurationException] shouldBe thrownBy {
      GoogleConfiguration(ConfigFactory.parseString(badKeyInUserMode))
    }

    val badKeyInServiceAccountMode =
      """
        |google {
        |  application-name = "cromwell"
        |
        |  auths = [
        |    {
        |      name = "name-service"
        |      scheme = "service_account"
        |      service-account-id-botched-key = "my-google-account"
        |      pem-file = "/yonder/file.pem"
        |    }
        |  ]
        |}
      """.stripMargin

    a[GoogleConfigurationException] shouldBe thrownBy {
      GoogleConfiguration(ConfigFactory.parseString(badKeyInServiceAccountMode))
    }
  }
}
