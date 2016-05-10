package cromwell.filesystems.gcs

import com.google.api.services.storage.StorageScopes
import com.typesafe.config.{Config, ConfigFactory}
import lenthall.config.ConfigValidationException
import lenthall.config.ValidatedConfig._
import org.slf4j.LoggerFactory

import scala.collection.JavaConverters._
import scala.language.postfixOps
import scalaz.Scalaz._
import scalaz.Validation.FlatMap._
import scalaz._


final case class GoogleConfiguration private (applicationName: String, authsByName: Map[String, GoogleAuthMode]) {
  def auth(name: String): ErrorOr[GoogleAuthMode] = {
    authsByName.get(name) match {
      case None =>
        val knownAuthNames = authsByName.keys.mkString(", ")
        s"`google` configuration stanza does not contain an auth named '$name'.  Known auth names: $knownAuthNames".failureNel
      case Some(a) => a.successNel
    }
  }
}

object GoogleConfiguration {

  private val log = LoggerFactory.getLogger("GoogleConfiguration")

  private val GoogleScopes = List(
    StorageScopes.DEVSTORAGE_FULL_CONTROL,
    StorageScopes.DEVSTORAGE_READ_WRITE,
    "https://www.googleapis.com/auth/genomics",
    "https://www.googleapis.com/auth/compute"
  )

  def apply(config: Config): GoogleConfiguration = {

    val googleConfig = config.getConfig("google")

    val appName = googleConfig.validateString("application-name")

    def buildAuth(authConfig: Config): ErrorOr[GoogleAuthMode] = {

      def serviceAccountAuth(authConfig: Config, name: String) = authConfig validateAny {
        cfg => ServiceAccountMode(name, cfg.getString("service-account-id"), cfg.getString("pem-file"), GoogleScopes)
      }

      def userAccountAuth(authConfig: Config, name: String) = authConfig validateAny {
        cfg => UserMode(name, cfg.getString("user"), cfg.getString("secrets-file"), cfg.getString("data-store-dir"), GoogleScopes)
      }

      def refreshTokenAuth(authConfig: Config, name: String) = authConfig validateAny {
        cfg => RefreshTokenMode(name, cfg.getString("client-id"), cfg.getString("client-secret"))
      }

      def applicationDefaultAuth(name: String) = ApplicationDefaultMode(name, GoogleScopes).successNel

      val name = authConfig.getString("name")
      val scheme = authConfig.getString("scheme")
      scheme match {
        case "service_account" => serviceAccountAuth(authConfig, name)
        case "user_account" => userAccountAuth(authConfig, name)
        case "refresh_token" => refreshTokenAuth(authConfig, name)
        case "application_default" => applicationDefaultAuth(name)
        case wut => s"Unsupported authentication scheme: $wut".failureNel
      }
    }

    val listOfErrorOrAuths: List[ErrorOr[GoogleAuthMode]] = googleConfig.getConfigList("auths").asScala.toList map buildAuth
    val errorOrAuthList: ErrorOr[List[GoogleAuthMode]] = listOfErrorOrAuths.sequence[ErrorOr, GoogleAuthMode]

    def uniqueAuthNames(list: List[GoogleAuthMode]): ErrorOr[Unit] = {
      val duplicateAuthNames = list.groupBy(_.name) collect { case (n, as) if as.size > 1 => n }
      if (duplicateAuthNames.nonEmpty) {
        ("Duplicate auth names: " + duplicateAuthNames.mkString(", ")).failureNel
      } else {
        ().successNel
      }
    }

    (appName |@| errorOrAuthList) { (_, _) } flatMap { case (name, list) =>
      uniqueAuthNames(list) map { _ =>
        GoogleConfiguration(name, list map { a => a.name -> a } toMap)
      }
    } match {
      case Success(r) => r
      case Failure(f) =>
        val errorMessages = f.list.mkString(", ")
        log.error(errorMessages)
        throw new ConfigValidationException("Google", errorMessages)
    }
  }
}
