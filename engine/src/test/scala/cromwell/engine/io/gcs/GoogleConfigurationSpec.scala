package cromwell.engine.io.gcs

import com.typesafe.config.{ConfigException, ConfigFactory}
import cromwell.util.ConfigUtil.ConfigValidationException
import org.scalatest.{Matchers, FlatSpec}

class GoogleConfigurationSpec extends FlatSpec with Matchers {

  behavior of "GoogleConfiguration"

  it should "parse a configuration stanza with service account" in {
    val serviceConfig =
      """
        |{
        | applicationName = "cromwell"
        |
        | cromwellAuthenticationScheme = "service_account"
        | serviceAuth = {
        |   serviceAccountId = "accountID"
        |   pemFile = "path/to/pem"
        | }
        |}
      """.stripMargin

    val gconf = GoogleConfiguration.fromConfig(ConfigFactory.parseString(serviceConfig)).get

    gconf.appName shouldBe "cromwell"
    gconf.cromwellAuthMode.isInstanceOf[ServiceAccountMode] shouldBe true
    val serviceAcccount = gconf.cromwellAuthMode.asInstanceOf[ServiceAccountMode]
    serviceAcccount.accountId shouldBe "accountID"
    serviceAcccount.pemPath shouldBe "path/to/pem"
  }

  it should "parse a configuration stanza with user account" in {
    val userConfig =
      """
        |{
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

    val gconf = GoogleConfiguration.fromConfig(ConfigFactory.parseString(userConfig)).get

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
        |{
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

    val gconf = GoogleConfiguration.fromConfig(ConfigFactory.parseString(refreshConfig)).get

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
        |{
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
       GoogleConfiguration.fromConfig(ConfigFactory.parseString(wrongConf)).get
    }
  }

  it should "not parse a configuration stanza with wrong cromwell auth" in {
    val unsupported =
      """
        |{
        | applicationName = "cromwell"
        |
        | cromwellAuthenticationScheme = "not supported"
        |}
      """.stripMargin

    a [ConfigValidationException] shouldBe thrownBy {
      GoogleConfiguration.fromConfig(ConfigFactory.parseString(unsupported)).get
    }

    val wrongSA =
      """
        |{
        | applicationName = "cromwell"
        |
        | cromwellAuthenticationScheme = "service_account"
        |}
      """.stripMargin

    a [ConfigException.Missing] shouldBe thrownBy {
      GoogleConfiguration.fromConfig(ConfigFactory.parseString(wrongSA)).get
    }

    val wrongSA2 =
      """
        |{
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
      GoogleConfiguration.fromConfig(ConfigFactory.parseString(wrongSA2)).get
    }
  }

  it should "not parse a configuration stanza with wrong user auth" in {
    val unsupported =
      """
        |{
        | applicationName = "cromwell"
        |
        | userAuthenticationScheme = "not supported"
        |}
      """.stripMargin

    a [ConfigValidationException] shouldBe thrownBy {
      GoogleConfiguration.fromConfig(ConfigFactory.parseString(unsupported)).get
    }

    val wrongUser =
      """
        |{
        | applicationName = "cromwell"
        |
        | cromwellAuthenticationScheme = "user_account"
        |}
      """.stripMargin

    a [ConfigException.Missing] shouldBe thrownBy {
      GoogleConfiguration.fromConfig(ConfigFactory.parseString(wrongUser)).get
    }

    val wrongUser2 =
      """
        |{
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
      GoogleConfiguration.fromConfig(ConfigFactory.parseString(wrongUser2)).get
    }
  }

}
