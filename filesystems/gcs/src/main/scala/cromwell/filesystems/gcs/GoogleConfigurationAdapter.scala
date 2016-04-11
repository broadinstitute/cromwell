package cromwell.filesystems.gcs

import com.typesafe.config.{Config, ConfigFactory}
import lenthall.config.ValidatedConfig._

import scala.util.Try
import scalaz.Scalaz._
import scalaz._

final case class GoogleConfigurationAdapter(appName: String, cromwellAuthMode: GoogleAuthMode, userAuthMode: Option[GoogleAuthMode]) {
  lazy val cromwellConf = {
    GoogleConfiguration(appName, cromwellAuthMode)
  }

  lazy val userConf = {
    userAuthMode map { c => GoogleConfiguration(appName, c) }
  }
}

object GoogleConfigurationAdapter {
  lazy val gcloudConf: Try[GoogleConfigurationAdapter] = Try(build())

  def build(conf: Config = ConfigFactory.load()): GoogleConfigurationAdapter = {
    val googleConf = conf.getConfig("google")

    val appName = googleConf.validateString("applicationName")

    def serviceAccountAuth = googleConf.getConfig("serviceAuth") validateAny {
       config => ServiceAccountMode(config.getString("serviceAccountId"), config.getString("pemFile"))
    }

    def userAccountAuth = googleConf.getConfig("userAuth") validateAny {
      config => UserMode(config.getString("user"), config.getString("secretsFile"), config.getString("dataStoreDir"))
    }

    val cromwellAuth = googleConf.validateString("cromwellAuthenticationScheme") match {
      case Success("service_account") => serviceAccountAuth
      case Success("user_account") => userAccountAuth
      case Success("application_default") => ApplicationDefaultMode.successNel
      case Success(unsupported) => s"Unsupported cromwellAuthenticationScheme: $unsupported".failureNel
      case Failure(f) => s"Could not find a value for cromwellAuthenticationScheme: $f".failureNel
    }

    def clientSecrets = googleConf.getConfig("refreshTokenAuth") validateAny {
      config => SimpleClientSecrets(config.getString("client_id"), config.getString("client_secret"))
    }

    val userAuth = googleConf.validateString("userAuthenticationScheme") match {
      case Success("refresh") => clientSecrets map { secrets => Option(RefreshTokenMode(secrets)) }
      case Success(unsupported) => s"Unsupported userAuthenticationMode: $unsupported".failureNel
      case Failure(_) => None.successNel
    }

    (appName |@| cromwellAuth |@| userAuth) {
      GoogleConfigurationAdapter(_, _, _)
    }  match {
      case Success(r) => r
      case Failure(f) => throw new Exception(s"""Invalid Google configuration: ${f.list.mkString("\n")}""")
    }
  }
}
