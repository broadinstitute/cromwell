package cromwell.filesystems.gcs

import cats.data.Validated._
import cats.instances.list._
import cats.syntax.cartesian._
import cats.syntax.traverse._
import cats.syntax.validated._
import com.google.api.services.storage.StorageScopes
import com.typesafe.config.Config
import cromwell.filesystems.gcs.auth._
import lenthall.config.ConfigValidationException
import lenthall.config.ValidatedConfig._
import cromwell.core.ErrorOr._
import org.slf4j.LoggerFactory

final case class GoogleConfiguration private (applicationName: String, authsByName: Map[String, GoogleAuthMode]) {

  def auth(name: String): ErrorOr[GoogleAuthMode] = {
    authsByName.get(name) match {
      case None =>
        val knownAuthNames = authsByName.keys.mkString(", ")
        s"`google` configuration stanza does not contain an auth named '$name'.  Known auth names: $knownAuthNames".invalidNel
      case Some(a) => a.validNel
    }
  }
}

object GoogleConfiguration {
  import scala.collection.JavaConverters._
  private val log = LoggerFactory.getLogger("GoogleConfiguration")

  val GoogleScopes = List(
    StorageScopes.DEVSTORAGE_FULL_CONTROL,
    StorageScopes.DEVSTORAGE_READ_WRITE,
    "https://www.googleapis.com/auth/genomics",
    "https://www.googleapis.com/auth/compute"
  ).asJava

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
        cfg => RefreshTokenMode(name, cfg.getString("client-id"), cfg.getString("client-secret"), GoogleScopes)
      }

      def applicationDefaultAuth(name: String): ErrorOr[GoogleAuthMode] = ApplicationDefaultMode(name).validNel

      val name = authConfig.getString("name")
      val scheme = authConfig.getString("scheme")
      scheme match {
        case "service_account" => serviceAccountAuth(authConfig, name)
        case "user_account" => userAccountAuth(authConfig, name)
        case "refresh_token" => refreshTokenAuth(authConfig, name)
        case "application_default" => applicationDefaultAuth(name)
        case wut => s"Unsupported authentication scheme: $wut".invalidNel
      }
    }

    val listOfErrorOrAuths: List[ErrorOr[GoogleAuthMode]] = googleConfig.getConfigList("auths").asScala.toList map buildAuth
    val errorOrAuthList: ErrorOr[List[GoogleAuthMode]] = listOfErrorOrAuths.sequence[ErrorOr, GoogleAuthMode]

    def uniqueAuthNames(list: List[GoogleAuthMode]): ErrorOr[Unit] = {
      val duplicateAuthNames = list.groupBy(_.name) collect { case (n, as) if as.size > 1 => n }
      if (duplicateAuthNames.nonEmpty) {
        ("Duplicate auth names: " + duplicateAuthNames.mkString(", ")).invalidNel
      } else {
        ().validNel
      }
    }

    (appName |@| errorOrAuthList) map { (_, _) } flatMap { case (name, list) =>
      uniqueAuthNames(list) map { _ =>
        GoogleConfiguration(name, list map { a => a.name -> a } toMap)
      }
    } match {
      case Valid(r) => r
      case Invalid(f) =>
        val errorMessages = f.toList.mkString(", ")
        log.error(errorMessages)
        throw new ConfigValidationException("Google", errorMessages)
    }
  }
}
