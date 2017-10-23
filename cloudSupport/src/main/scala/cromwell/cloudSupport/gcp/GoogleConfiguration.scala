package cromwell.cloudsupport.gcp

import java.io.IOException

import cats.data.Validated._
import cats.instances.list._
import cats.syntax.traverse._
import cats.syntax.validated._
import com.google.api.client.googleapis.auth.oauth2.GoogleCredential
import com.google.api.client.http.{HttpRequest, HttpRequestInitializer}
import com.google.api.services.storage.StorageScopes
import com.typesafe.config.{Config, ConfigException}
import cromwell.cloudsupport.gcp.auth.ServiceAccountMode.{JsonFileFormat, PemFileFormat}
import cromwell.cloudsupport.gcp.auth._
import common.exception.MessageAggregation
import common.validation.ErrorOr._
import common.validation.Validation._
import net.ceedubs.ficus.Ficus._
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
  import scala.concurrent.duration._
  import scala.language.postfixOps

  lazy val DefaultConnectionTimeout = 3 minutes
  lazy val DefaultReadTimeout = 3 minutes

  lazy val DefaultRequestInitializer = GoogleConfiguration.withCustomTimeouts(new GoogleCredential.Builder().build())

  def withCustomTimeouts(httpRequestInitializer: HttpRequestInitializer,
                         connectionTimeout: FiniteDuration = DefaultConnectionTimeout,
                         readTimeout: FiniteDuration = DefaultReadTimeout) = {
    new HttpRequestInitializer() {
      @throws[IOException]
      override def initialize(httpRequest: HttpRequest) = {
        httpRequestInitializer.initialize(httpRequest)
        httpRequest.setConnectTimeout(connectionTimeout.toMillis.toInt)
        httpRequest.setReadTimeout(readTimeout.toMillis.toInt)
        ()
      }
    }
  }

  private val log = LoggerFactory.getLogger("GoogleConfiguration")

  final case class GoogleConfigurationException(errorMessages: List[String]) extends MessageAggregation {
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
        (authConfig.getAs[String]("pem-file"), authConfig.getAs[String]("json-file")) match {
          case (Some(pem), None) => ServiceAccountMode(name, PemFileFormat(authConfig.as[String]("service-account-id"), pem), GoogleScopes)
          case (None, Some(json)) => ServiceAccountMode(name, JsonFileFormat(json), GoogleScopes)
          case (None, None) => throw new ConfigException.Generic(s"""No credential configuration was found for service account "$name". See reference.conf under the google.auth, service-account section for supported credential formats.""")
          case (Some(_), Some(_)) => throw new ConfigException.Generic(s"""Both a pem file and a json file were supplied for service account "$name" in the configuration file. Only one credential file can be supplied for the same service account. Please choose between the two.""")
        }
      }

      def userAccountAuth(authConfig: Config, name: String): ErrorOr[GoogleAuthMode] =  validate {
        UserMode(name, authConfig.as[String]("user"), authConfig.as[String]("secrets-file"), authConfig.as[String]("data-store-dir"), GoogleScopes)
      }

      def refreshTokenAuth(authConfig: Config, name: String): ErrorOr[GoogleAuthMode] = validate {
        RefreshTokenMode(name, authConfig.as[String]("client-id"), authConfig.as[String]("client-secret"), GoogleScopes)
      }

      def applicationDefaultAuth(name: String): ErrorOr[GoogleAuthMode] = ApplicationDefaultMode(name).validNel

      def userServiceAccountAuth(name: String): ErrorOr[GoogleAuthMode] = validate {
        UserServiceAccountMode(name, GoogleScopes)
      }

      val name = authConfig.getString("name")
      val scheme = authConfig.getString("scheme")
      scheme match {
        case "service_account" => serviceAccountAuth(authConfig, name)
        case "user_account" => userAccountAuth(authConfig, name)
        case "refresh_token" => refreshTokenAuth(authConfig, name)
        case "application_default" => applicationDefaultAuth(name)
        case "user_service_account" => userServiceAccountAuth(name)
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

    (appName, errorOrAuthList) flatMapN { (name, list) =>
      uniqueAuthNames(list) map { _ =>
        GoogleConfiguration(name, list map { a => a.name -> a } toMap)
      }
    } match {
      case Valid(r) => r
      case Invalid(f) =>
        val errorMessages = f.toList.mkString(", ")
        log.error(errorMessages)
        throw GoogleConfigurationException(f.toList)
    }
  }
}
