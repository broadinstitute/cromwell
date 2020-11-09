package cromwell.cloudsupport.gcp

import java.net.URL

import better.files.File
import cats.syntax.all._
import com.google.api.client.http.GenericUrl
import com.google.api.client.testing.http.MockHttpTransport
import com.typesafe.config.{ConfigException, ConfigFactory}
import common.assertion.CromwellTimeoutSpec
import cromwell.cloudsupport.gcp.GoogleConfiguration.GoogleConfigurationException
import cromwell.cloudsupport.gcp.auth.ServiceAccountMode.{JsonFileFormat, PemFileFormat}
import cromwell.cloudsupport.gcp.auth._
import org.scalatest.flatspec.AnyFlatSpec
import org.scalatest.matchers.should.Matchers


class GoogleConfigurationSpec extends AnyFlatSpec with CromwellTimeoutSpec with Matchers {

  behavior of "GoogleConfiguration"

  it should "parse all manner of well-formed auths" in {
    val pemMockFile = File.newTemporaryFile()
    val jsonMockFile = File.newTemporaryFile()

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
        |      secrets-file = "${pemMockFile.pathAsString}"
        |      data-store-dir = "/where/the/data/at"
        |    },
        |    {
        |      name = "name-pem-service"
        |      scheme = "service_account"
        |      service-account-id = "my-google-account"
        |      pem-file = "${pemMockFile.pathAsString}"
        |    },
        |    {
        |      name = "name-json-service"
        |      scheme = "service_account"
        |      json-file = "${jsonMockFile.pathAsString}"
        |    },
        |    {
        |      name = "name-user-service-account"
        |      scheme = "user_service_account"
        |    }
        |  ]
        |}
        |
      """.stripMargin

    val gconf = GoogleConfiguration(ConfigFactory.parseString(righteousGoogleConfig))

    gconf.applicationName shouldBe "cromwell"
    gconf.authsByName should have size 6

    val auths = gconf.authsByName.values

    val appDefault = (auths collectFirst { case a: ApplicationDefaultMode => a }).get
    appDefault.name shouldBe "name-default"

    val refreshToken = (auths collectFirst { case a: RefreshTokenMode => a }).get
    refreshToken.name shouldBe "name-refresh"
    refreshToken.clientSecret shouldBe "secret_secret"
    refreshToken.clientId shouldBe "secret_id"

    val userServiceAccount = (auths collectFirst { case a: UserServiceAccountMode => a }).get
    userServiceAccount.name shouldBe "name-user-service-account"

    val user = (auths collectFirst { case a: UserMode => a }).get
    user.name shouldBe "name-user"
    user.secretsPath shouldBe pemMockFile.pathAsString

    val servicePem = (auths collectFirst { case a: ServiceAccountMode if a.name == "name-pem-service" => a }).get
    servicePem.name shouldBe "name-pem-service"
    servicePem.fileFormat.asInstanceOf[PemFileFormat].accountId shouldBe "my-google-account"
    servicePem.fileFormat.file shouldBe pemMockFile.pathAsString

    val serviceJson = (auths collectFirst { case a: ServiceAccountMode if a.name == "name-json-service" => a }).get
    serviceJson.name shouldBe "name-json-service"
    serviceJson.fileFormat.isInstanceOf[JsonFileFormat] shouldBe true
    serviceJson.fileFormat.file shouldBe jsonMockFile.pathAsString

    pemMockFile.delete(swallowIOExceptions = true)
    jsonMockFile.delete(swallowIOExceptions = true)
  }

  it should "return a known auth" in {
    val config =
      """|google {
         |  application-name = "cromwell"
         |
         |  auths = [
         |    {
         |      name = "name-default"
         |      scheme = "application_default"
         |    }
         |  ]
         |}
         |""".stripMargin

    val googleConfiguration = GoogleConfiguration(ConfigFactory.parseString(config))
    googleConfiguration.auth("name-default").map(_.name) should be("name-default".valid)
  }

  it should "not return an unknown auth" in {
    val config =
      """|google {
         |  application-name = "cromwell"
         |
         |  auths = [
         |    {
         |      name = "name-default"
         |      scheme = "application_default"
         |    }
         |  ]
         |}
         |""".stripMargin

    val googleConfiguration = GoogleConfiguration(ConfigFactory.parseString(config))
    googleConfiguration.auth("name-botched") should be(
      "`google` configuration stanza does not contain an auth named 'name-botched'.  Known auth names: name-default"
        .invalidNel)
  }

  it should "create an initializer with custom timeouts" in {
    val transport = new MockHttpTransport()
    val initializer = GoogleConfiguration.withCustomTimeouts(request => {
      request.getHeaders.set("custom_init", "ok")
      ()
    })
    val factory = transport.createRequestFactory(initializer)
    val request = factory.buildGetRequest(new GenericUrl(new URL("http://example.com")))
    request.getConnectTimeout should be(180000)
    request.getReadTimeout should be(180000)
    request.getHeaders.get("custom_init") should be("ok")
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

    the[GoogleConfigurationException] thrownBy {
      GoogleConfiguration(ConfigFactory.parseString(applessGoogleConfig))
    } should have message "Google configuration:\nString: 2: No configuration setting found for key 'application-name'"
  }

  it should "not parse a configuration stanza with double service account credentials" in {
    val doubleServiceAccountCredentials =
      """
        |google {
        |  application-name = "cromwell"
        |
        |  auths = [
        |    {
        |      name = "service-account"
        |      scheme = "service_account"
        |      service-account-id = "my-google-account"
        |      pem-file = "path/to/file.pem"
        |      json-file = "path/to/json.pem"
        |    }
        |  ]
        |}
      """.stripMargin

    the[GoogleConfigurationException] thrownBy {
      GoogleConfiguration(ConfigFactory.parseString(doubleServiceAccountCredentials))
    } should have message "Google configuration:\n" +
      "Both a pem file and a json file were supplied for service account \"service-account\" in the configuration " +
      "file. Only one credential file can be supplied for the same service account. Please choose between the two."
  }

  it should "not parse a configuration stanza without service account credentials" in {
    val noServiceAccountCredentials =
      """
        |google {
        |  application-name = "cromwell"
        |
        |  auths = [
        |    {
        |      name = "service-account"
        |      scheme = "service_account"
        |      service-account-id = "my-google-account"
        |    }
        |  ]
        |}
      """.stripMargin

    the[GoogleConfigurationException] thrownBy {
      GoogleConfiguration(ConfigFactory.parseString(noServiceAccountCredentials))
    } should have message "Google configuration:\n" +
      "No credential configuration was found for service account \"service-account\". See reference.conf under the " +
      "google.auth, service-account section for supported credential formats."
  }

  it should "not parse a configuration stanza with an unsupported authentication scheme" in {
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

    the[GoogleConfigurationException] thrownBy {
      GoogleConfiguration(ConfigFactory.parseString(unsupported))
    } should have message "Google configuration:\nUnsupported authentication scheme: not supported"
  }

  it should "not parse a configuration stanza without a schema" in {
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

    the[ConfigException.Missing] thrownBy {
      GoogleConfiguration(ConfigFactory.parseString(schemeless))
    } should have message "String: 6: No configuration setting found for key 'scheme'"
  }

  it should "not parse a configuration stanza without an auth name" in {
    val nameless =
      """
        |google {
        |  application-name = "cromwell"
        |
        |  auths = [
        |    {
        |      scheme = "application_default"
        |    }
        |  ]
        |}
      """.stripMargin

    the[ConfigException.Missing] thrownBy {
      GoogleConfiguration(ConfigFactory.parseString(nameless))
    } should have message "String: 6: No configuration setting found for key 'name'"
  }

  it should "not parse a configuration stanza with a bad client-id in refresh token mode" in {
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

    the[GoogleConfigurationException] thrownBy {
      GoogleConfiguration(ConfigFactory.parseString(badKeyInRefreshTokenMode))
    } should have message "Google configuration:\nString: 6: No configuration setting found for key 'client-id'"
  }

  it should "parse a configuration stanza without a user in user mode" in {
    val config =
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

    val googleConfiguration = GoogleConfiguration(ConfigFactory.parseString(config))
    googleConfiguration.auth("name-user").map(_.name) should be("name-user".valid)
  }

  it should "not parse a configuration stanza without a service-account-id in service account mode" in {
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

    the[GoogleConfigurationException] thrownBy {
      GoogleConfiguration(ConfigFactory.parseString(badKeyInServiceAccountMode))
    } should have message "Google configuration:\nString: 6: No configuration setting found for key 'service-account-id'"
  }

  it should "not parse a configuration stanza with a duplicate auth name" in {
    val duplicateAuthName =
      """|google {
         |  application-name = "cromwell"
         |
         |  auths = [
         |    {
         |      name = "name-default"
         |      scheme = "application_default"
         |    }
         |    {
         |      name = "name-default"
         |      scheme = "application_default"
         |    }
         |  ]
         |}
         |""".stripMargin

    the[GoogleConfigurationException] thrownBy {
      GoogleConfiguration(ConfigFactory.parseString(duplicateAuthName))
    } should have message "Google configuration:\nDuplicate auth names: name-default"
  }
}
