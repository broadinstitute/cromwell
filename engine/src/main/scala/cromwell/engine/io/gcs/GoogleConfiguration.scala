package cromwell.engine.io.gcs

import com.typesafe.config.Config
import cromwell.util.ConfigUtil._
import org.slf4j.LoggerFactory

import scala.language.postfixOps
import scala.util.Try
import scalaz.Scalaz._
import scalaz._

// Google Authentication modes supported for Cromwell
sealed trait GoogleCromwellAuthMode
final case class ServiceAccountMode(accountId: String, pemPath: String) extends GoogleCromwellAuthMode
final case class UserMode(user: String, secretsFile: String, datastoreDir: String) extends GoogleCromwellAuthMode
case object ApplicationDefaultMode extends GoogleCromwellAuthMode

// Google Authentication modes supported for User
sealed trait GoogleUserAuthMode
final case class Refresh(clientSecrets: ClientSecrets) extends GoogleUserAuthMode

trait ClientSecrets {
  val clientId: String
  val clientSecret: String
}

final case class SimpleClientSecrets(clientId: String, clientSecret: String) extends ClientSecrets
final case class GoogleConfiguration(appName: String, cromwellAuthMode: GoogleCromwellAuthMode, userAuthMode: Option[GoogleUserAuthMode])

object GoogleConfiguration {

  private val log = LoggerFactory.getLogger("GoogleConfiguration")

  def fromConfig(config: Config): Try[GoogleConfiguration] = Try {

    val appName = config.validateString("applicationName")

    def serviceAccountAuth = config.getConfig("serviceAuth") validateAny {
      cfg => ServiceAccountMode(cfg.getString("serviceAccountId"), cfg.getString("pemFile"))
    }

    def userAccountAuth = config.getConfig("userAuth") validateAny {
      cfg => UserMode(cfg.getString("user"), cfg.getString("secretsFile"), cfg.getString("dataStoreDir"))
    }

    val cromwellAuth = config.validateString("cromwellAuthenticationScheme") match {
      case Success("service_account") => serviceAccountAuth
      case Success("user_account") => userAccountAuth
      case Success("application_default") => ApplicationDefaultMode.successNel
      case Success(unsupported) => s"Unsupported cromwellAuthenticationScheme: $unsupported".failureNel
      case Failure(f) => s"Could not find a value for cromwellAuthenticationScheme: $f".failureNel
    }

    def clientSecrets = config.getConfig("refreshTokenAuth") validateAny {
      cfg => SimpleClientSecrets(cfg.getString("client_id"), cfg.getString("client_secret"))
    }

    val userAuth = config.validateString("userAuthenticationScheme") match {
      case Success("refresh") => clientSecrets map { secrets => Option(Refresh(secrets)) }
      case Success(unsupported) => s"Unsupported userAuthenticationMode: $unsupported".failureNel
      case Failure(_) => None.successNel
    }

    (appName |@| cromwellAuth |@| userAuth) {
      GoogleConfiguration(_, _, _)
    }  match {
      case Success(r) => r
      case Failure(f) =>
        val errorMessages = f.list.mkString(", ")
        log.error(errorMessages)
        throw new ConfigValidationException("Google", errorMessages)
    }
  }
}
