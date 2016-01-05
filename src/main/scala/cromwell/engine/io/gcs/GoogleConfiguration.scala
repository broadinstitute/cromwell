package cromwell.engine.io.gcs

import com.typesafe.config.{Config, ConfigFactory}
import cromwell.util.ConfigUtil._

import scala.util.Try
import scalaz.Scalaz._
import scalaz._


// Google Authentication modes supported for Cromwell
sealed trait GoogleCromwellAuthMode
final case class ServiceAccountMode(accountId: String, pemPath: String) extends GoogleCromwellAuthMode
final case class UserMode(user: String, secretsFile: String, datastoreDir: String) extends GoogleCromwellAuthMode

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

  lazy val gcloudConf: Try[GoogleConfiguration] = Try(build())

  def build(conf: Config = ConfigFactory.load()): GoogleConfiguration = {
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
      case Success(unsupported) => s"Unsupported cromwellAuthenticationScheme: $unsupported".failureNel
      case Failure(f) => s"Could not find a value for cromwellAuthenticationScheme: $f".failureNel
    }

    def clientSecrets = googleConf.getConfig("refreshTokenAuth") validateAny {
      config => SimpleClientSecrets(config.getString("client_id"), config.getString("client_secret"))
    }
    val userAuth = googleConf.validateString("userAuthenticationScheme") match {
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
        throw new ConfigValidationException("Google", errorMessages)
    }
  }
}
