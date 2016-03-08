package cromwell.engine.io.gcs

import com.typesafe.config.{ConfigException, ConfigFactory}
import cromwell.util.ConfigUtil.ConfigValidationException
import org.scalatest.{Matchers, FlatSpec}

class GoogleConfigurationSpec extends FlatSpec with Matchers {

  behavior of "GoogleConfiguration"

  it should "parse a configuration stanza with service account" in {
    val serviceConfig =
      """
        |google {
        | applicationName = "cromwell"
        |
        | cromwellAuthenticationScheme = "service_account"
        | serviceAuth = {
        |   serviceAccountId = "accountID"
        |   pemFile = "path/to/pem"
        | }
        |}
      """.stripMargin

    val gconf = GoogleConfiguration.build(ConfigFactory.parseString(serviceConfig))

    gconf.appName shouldBe "cromwell"
    gconf.cromwellAuthMode.isInstanceOf[ServiceAccountMode] shouldBe true
    val serviceAcccount = gconf.cromwellAuthMode.asInstanceOf[ServiceAccountMode]
    serviceAcccount.accountId shouldBe "accountID"
    serviceAcccount.pemPath shouldBe "path/to/pem"
  }

  it should "parse a configuration stanza with user account" in {
    val userConfig =
      """
        |google {
        | applicationName = "cromwell"
        |
        | cromwellAuthenticationScheme = "user_account"
        | userAuth {
        |    user = "username"
        |    secretsFile = "path/to/secrets"
        |    dataStoreDir = "path/to/store"
        |  }
        |}
      """.stripMargin

    val gconf = GoogleConfiguration.build(ConfigFactory.parseString(userConfig))

    gconf.appName shouldBe "cromwell"
    gconf.cromwellAuthMode.isInstanceOf[UserMode] shouldBe true
    val userAccount = gconf.cromwellAuthMode.asInstanceOf[UserMode]
    userAccount.user shouldBe "username"
    userAccount.secretsFile shouldBe "path/to/secrets"
    userAccount.datastoreDir shouldBe "path/to/store"
  }

  it should "parse a configuration stanza with refresh token scheme" in {
    val refreshConfig =
      """
        |google {
        | applicationName = "cromwell"
        |
        | cromwellAuthenticationScheme = "service_account"
        | serviceAuth = {
        |   serviceAccountId = "accountID"
        |   pemFile = "path/to/pem"
        | }
        |
        |  userAuthenticationScheme = "refresh"
        |  refreshTokenAuth = {
        |    client_id = "clientID"
        |    client_secret = "clientSecret"
        |  }
        |}
      """.stripMargin

    val gconf = GoogleConfiguration.build(ConfigFactory.parseString(refreshConfig))

    val refreshConf = gconf.userAuthMode
    refreshConf shouldBe defined
    refreshConf.get.isInstanceOf[Refresh] shouldBe true
    val refresh = refreshConf.get.asInstanceOf[Refresh]
    refresh.clientSecrets.clientId shouldBe "clientID"
    refresh.clientSecrets.clientSecret shouldBe "clientSecret"
  }

  it should "not parse a configuration stanza without applicationName" in {
    val wrongConf =
      """
        |google {
        |
        | cromwellAuthenticationScheme = "user_account"
        | userAuth {
        |    user = "username"
        |    secretsFile = "path/to/secrets"
        |    dataStoreDir = "path/to/store"
        |  }
        |}
      """.stripMargin

    a [ConfigValidationException] shouldBe thrownBy {
       GoogleConfiguration.build(ConfigFactory.parseString(wrongConf))
    }
  }

  it should "not parse a configuration stanza with wrong cromwell auth" in {
    val unsupported =
      """
        |google {
        | applicationName = "cromwell"
        |
        | cromwellAuthenticationScheme = "not supported"
        |}
      """.stripMargin

    a [ConfigValidationException] shouldBe thrownBy {
      GoogleConfiguration.build(ConfigFactory.parseString(unsupported))
    }

    val wrongSA =
      """
        |google {
        | applicationName = "cromwell"
        |
        | cromwellAuthenticationScheme = "service_account"
        |}
      """.stripMargin

    a [ConfigException.Missing] shouldBe thrownBy {
      GoogleConfiguration.build(ConfigFactory.parseString(wrongSA))
    }

    val wrongSA2 =
      """
        |google {
        | applicationName = "cromwell"
        |
        | cromwellAuthenticationScheme = "service_account"
        | serviceAuth = {
        |   wrongKey = "accountID"
        |   pemFile = "path/to/pem"
        | }
        |}
      """.stripMargin

    a [ConfigException.Missing] shouldBe thrownBy {
      GoogleConfiguration.build(ConfigFactory.parseString(wrongSA2))
    }
  }

  it should "not parse a configuration stanza with wrong user auth" in {
    val unsupported =
      """
        |google {
        | applicationName = "cromwell"
        |
        | userAuthenticationScheme = "not supported"
        |}
      """.stripMargin

    a [ConfigValidationException] shouldBe thrownBy {
      GoogleConfiguration.build(ConfigFactory.parseString(unsupported))
    }

    val wrongUser =
      """
        |google {
        | applicationName = "cromwell"
        |
        | cromwellAuthenticationScheme = "user_account"
        |}
      """.stripMargin

    a [ConfigException.Missing] shouldBe thrownBy {
      GoogleConfiguration.build(ConfigFactory.parseString(wrongUser))
    }

    val wrongUser2 =
      """
        |google {
        | applicationName = "cromwell"
        |
        |cromwellAuthenticationScheme = "user_account"
        | userAuth {
        |    wrongKey = "username"
        |    secretsFile = "path/to/secrets"
        |    dataStoreDir = "path/to/store"
        |  }
        |}
      """.stripMargin

    a [ConfigException.Missing] shouldBe thrownBy {
      GoogleConfiguration.build(ConfigFactory.parseString(wrongUser2))
    }
  }

}
