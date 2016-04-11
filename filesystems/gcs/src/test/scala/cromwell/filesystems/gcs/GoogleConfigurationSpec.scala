package cromwell.filesystems.gcs

import com.typesafe.config.{ConfigException, ConfigFactory}
import lenthall.config.ConfigValidationException
import org.scalatest.{FlatSpec, Matchers}

class GoogleConfigurationSpec extends FlatSpec with Matchers {

  behavior of "GoogleConfiguration"

  it should "parse a configuration stanza with service account" in {
    val serviceConfig =
      """
        |applicationName = "cromwell"
        |
        |authenticationScheme = "service_account"
        |serviceAuth = {
        |  serviceAccountId = "accountID"
        |  pemFile = "path/to/pem"
        |}
      """.stripMargin

    val gconf = GoogleConfiguration(ConfigFactory.parseString(serviceConfig))

    gconf.appName shouldBe "cromwell"
    gconf.authMode.isInstanceOf[ServiceAccountMode] shouldBe true
    val serviceAcccount = gconf.authMode.asInstanceOf[ServiceAccountMode]
    serviceAcccount.accountId shouldBe "accountID"
    serviceAcccount.pemPath shouldBe "path/to/pem"
  }

  it should "parse a configuration stanza with user account" in {
    val userConfig =
      """
        |applicationName = "cromwell"
        |
        |authenticationScheme = "user_account"
        |userAuth {
        |   user = "username"
        |   secretsFile = "path/to/secrets"
        |   dataStoreDir = "path/to/store"
        | }
      """.stripMargin

    val gconf = GoogleConfiguration(ConfigFactory.parseString(userConfig))

    gconf.appName shouldBe "cromwell"
    gconf.authMode.isInstanceOf[UserMode] shouldBe true
    val userAccount = gconf.authMode.asInstanceOf[UserMode]
    userAccount.user shouldBe "username"
    userAccount.secretsFile shouldBe "path/to/secrets"
    userAccount.datastoreDir shouldBe "path/to/store"
  }

  it should "parse a configuration stanza with refresh token scheme" in {
    val refreshConfig =
      """
        |applicationName = "cromwell"

        | authenticationScheme = "refresh_token"
        | refreshTokenAuth = {
        |   client_id = "clientID"
        |   client_secret = "clientSecret"
        | }
      """.stripMargin

    val gconf = GoogleConfiguration(ConfigFactory.parseString(refreshConfig))

    val refreshConf = gconf.authMode
    gconf.authMode.isInstanceOf[RefreshTokenMode] shouldBe true
    val refresh = gconf.authMode.asInstanceOf[RefreshTokenMode]
    refresh.clientSecrets.clientId shouldBe "clientID"
    refresh.clientSecrets.clientSecret shouldBe "clientSecret"
  }

  it should "not parse a configuration stanza without applicationName" in {
    val wrongConf =
      """
        |authenticationScheme = "user_account"
        |userAuth {
        |   user = "username"
        |   secretsFile = "path/to/secrets"
        |   dataStoreDir = "path/to/store"
        | }
      """.stripMargin

    a [ConfigValidationException] shouldBe thrownBy {
       GoogleConfiguration(ConfigFactory.parseString(wrongConf))
    }
  }

  it should "not parse a configuration stanza with wrong cromwell auth" in {
    val unsupported =
      """
        |applicationName = "cromwell"
        |
        |authenticationScheme = "not supported"
      """.stripMargin

    a [ConfigValidationException] shouldBe thrownBy {
      GoogleConfiguration(ConfigFactory.parseString(unsupported))
    }

    val wrongSA =
      """
        |applicationName = "cromwell"
        |
        |authenticationScheme = "service_account"
      """.stripMargin

    a [ConfigException.Missing] shouldBe thrownBy {
      GoogleConfiguration(ConfigFactory.parseString(wrongSA))
    }

    val wrongSA2 =
      """
        |applicationName = "cromwell"
        |
        |authenticationScheme = "service_account"
        |serviceAuth = {
        |  wrongKey = "accountID"
        |  pemFile = "path/to/pem"
        |}
      """.stripMargin

    a [ConfigException.Missing] shouldBe thrownBy {
      GoogleConfiguration(ConfigFactory.parseString(wrongSA2))
    }
  }

  it should "not parse a configuration stanza with wrong user auth" in {

    val wrongUser =
      """
        |applicationName = "cromwell"
        |
        |authenticationScheme = "user_account"
      """.stripMargin

    a [ConfigException.Missing] shouldBe thrownBy {
      GoogleConfiguration(ConfigFactory.parseString(wrongUser))
    }

    val wrongUser2 =
      """
        |applicationName = "cromwell"
        |
        |authenticationScheme = "user_account"
        |userAuth {
        |   wrongKey = "username"
        |   secretsFile = "path/to/secrets"
        |   dataStoreDir = "path/to/store"
        |}
      """.stripMargin

    a [ConfigException.Missing] shouldBe thrownBy {
      GoogleConfiguration(ConfigFactory.parseString(wrongUser2))
    }
  }

}
