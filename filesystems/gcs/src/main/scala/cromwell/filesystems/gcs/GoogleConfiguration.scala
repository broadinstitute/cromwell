package cromwell.filesystems.gcs

import cats.data.Validated._
import cats.instances.list._
import cats.syntax.cartesian._
import cats.syntax.traverse._
import cats.syntax.validated._
import com.google.api.services.storage.StorageScopes
import com.typesafe.config.Config
import cromwell.filesystems.gcs.auth._
import lenthall.exception.MessageAggregation
import lenthall.validation.ErrorOr._
import lenthall.validation.Validation._
import org.slf4j.LoggerFactory
import net.ceedubs.ficus.Ficus._

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

  case class GoogleConfigurationException(errorMessages: List[String]) extends MessageAggregation {
    override val exceptionContext = "Google configuration"
  }

  val GoogleScopes = List(
    StorageScopes.DEVSTORAGE_FULL_CONTROL,
    StorageScopes.DEVSTORAGE_READ_WRITE,
    "https://www.googleapis.com/auth/genomics",
    "https://www.googleapis.com/auth/compute"
  ).asJava
  
  def apply(config: Config): GoogleConfiguration = {

    val googleConfig = config.getConfig("google")

    val appName = validate { googleConfig.as[String]("application-name") }

    def buildAuth(authConfig: Config): ErrorOr[GoogleAuthMode] = {

      def serviceAccountAuth(authConfig: Config, name: String): ErrorOr[GoogleAuthMode] = validate {
        ServiceAccountMode(name, authConfig.as[String]("service-account-id"), authConfig.as[String]("pem-file"), GoogleScopes)
      }

      def userAccountAuth(authConfig: Config, name: String): ErrorOr[GoogleAuthMode] =  validate {
        UserMode(name, authConfig.as[String]("user"), authConfig.as[String]("secrets-file"), authConfig.as[String]("data-store-dir"), GoogleScopes)
      }

      def refreshTokenAuth(authConfig: Config, name: String): ErrorOr[GoogleAuthMode] = validate {
        RefreshTokenMode(name, authConfig.as[String]("client-id"), authConfig.as[String]("client-secret"), GoogleScopes)
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

    val listOfErrorOrAuths: List[ErrorOr[GoogleAuthMode]] = googleConfig.as[List[Config]]("auths") map buildAuth
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
        throw new GoogleConfigurationException(f.toList)
    }
  }
}
