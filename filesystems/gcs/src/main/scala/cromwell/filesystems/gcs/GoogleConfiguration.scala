package cromwell.filesystems.gcs

import com.typesafe.config.Config
import lenthall.config.ConfigValidationException
import lenthall.config.ValidatedConfig._

import scalaz.Scalaz._
import scalaz._

// Google Authentication modes supported for Cromwell
//FIXME not sealed yet because used in GoogleConfigurationAdapter temporarily
trait GoogleAuthMode
final case class ServiceAccountMode(accountId: String, pemPath: String) extends GoogleAuthMode
final case class UserMode(user: String, secretsFile: String, datastoreDir: String) extends GoogleAuthMode
final case class RefreshTokenMode(clientSecrets: ClientSecrets) extends GoogleAuthMode
case object ApplicationDefaultMode extends GoogleAuthMode

trait ClientSecrets {
  val clientId: String
  val clientSecret: String
}


final case class GoogleConfiguration(appName: String, authMode: GoogleAuthMode)

object GoogleConfiguration {
  val RefreshTokenOptionKey = "refresh_token"

  def apply(googleConf: Config): GoogleConfiguration = {

    val appName = googleConf.validateString("applicationName")

    def serviceAccountAuth = googleConf.getConfig("serviceAuth") validateAny {
      config => ServiceAccountMode(config.getString("serviceAccountId"), config.getString("pemFile"))
    }

    def userAccountAuth = googleConf.getConfig("userAuth") validateAny {
      config => UserMode(config.getString("user"), config.getString("secretsFile"), config.getString("dataStoreDir"))
    }

    def refreshTokenAuth = googleConf.getConfig("refreshTokenAuth") validateAny {
      config => RefreshTokenMode(SimpleClientSecrets(config.getString("client_id"), config.getString("client_secret")))
    }

    val cromwellAuth = googleConf.validateString("authenticationScheme") match {
      case Success("service_account") => serviceAccountAuth
      case Success("user_account") => userAccountAuth
      case Success("refresh_token") => refreshTokenAuth
      case Success("application_default") => ApplicationDefaultMode.successNel
      case Success(unsupported) => s"Unsupported cromwellAuthenticationScheme: $unsupported".failureNel
      case Failure(f) => f.failure
    }

    (appName |@| cromwellAuth) {
      GoogleConfiguration(_, _)
    }  match {
      case Success(r) => r
      case Failure(f) => throw new ConfigValidationException("Google", f)
    }
  }
}
