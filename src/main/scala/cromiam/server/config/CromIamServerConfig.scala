package cromiam.server.config

import akka.http.scaladsl.settings.ServerSettings
import cats.syntax.validated._
import cats.syntax.cartesian._
import com.typesafe.config.{Config, ConfigFactory}
import lenthall.validation.ErrorOr.ErrorOr

import scala.collection.JavaConverters._
import scala.util.{Failure, Success, Try}
import cromiam.server.config.CromIamServerConfig._

final case class CromIamServerConfig(cromIamConfig: CromIamConfig,
                                     cromwellConfig: ServiceConfig,
                                     samConfig: ServiceConfig,
                                     swaggerOauthConfig: SwaggerOauthConfig)

object CromIamServerConfig {
  def getFromConfig(conf: Config): ErrorOr[CromIamServerConfig] = {
    val cromIamConfig = CromIamConfig.getFromConfig(conf, "cromiam")
    val cromwellConfig = ServiceConfig.getFromConfig(conf, "cromwell")
    val samConfig = ServiceConfig.getFromConfig(conf, "sam")
    val googleConfig = SwaggerOauthConfig.getFromConfig(conf, "swagger_oauth")

    (cromIamConfig |@| cromwellConfig |@| samConfig |@| googleConfig) map CromIamServerConfig.apply
  }

  private[config] def getValidatedConfigPath[A](typename: String, conf: Config, path: String, getter: (Config, String) => A, default: Option[A] = None): ErrorOr[A] = {
    if (conf.hasPath(path)) {
      Try(getter.apply(conf, path)) match {
        case Success(s) => s.validNel
        case Failure(e) => s"Unable to read valid value at '$path': ${e.getMessage}".invalidNel
      }
    } else default match {
      case Some(d) => d.validNel
      case None => s"Configuration does not have path $path".invalidNel
    }
  }

  private[config] implicit final class ValidatingConfig(val conf: Config) extends AnyVal {
    def getValidatedString(path: String, default: Option[String] = None): ErrorOr[String] = getValidatedConfigPath("string", conf, path, (c, p) => c.getString(p), default)
    def getValidatedInt(path: String): ErrorOr[Int] = getValidatedConfigPath("integer", conf, path, (c, p) => c.getInt(p))
    def getValidatedStringList(path: String): ErrorOr[List[String]] = getValidatedConfigPath[List[String]]("string list", conf, path, (c, p) => c.getStringList(p).asScala.toList)
  }
}

final case class CromIamConfig(http: ServiceConfig, serverSettings: ServerSettings)

object CromIamConfig {

  private def getValidatedServerSettings: ErrorOr[ServerSettings] = Try(ServerSettings(ConfigFactory.load())) match {
    case Success(serverSettings) => serverSettings.validNel
    case Failure(e) => s"Unable to generate server settings from configuration file: ${e.getMessage}".invalidNel
  }

  private[config] def getFromConfig(conf: Config, basePath: String): ErrorOr[CromIamConfig] = {
    val serviceConfig = ServiceConfig.getFromConfig(conf, s"$basePath")
    val serverSettings = getValidatedServerSettings

    (serviceConfig |@| serverSettings) map CromIamConfig.apply
  }
}

final case class ServiceConfig(interface: String, port: Int, scheme: String)

object ServiceConfig {
  private[config] def getFromConfig(conf: Config, basePath: String): ErrorOr[ServiceConfig] = {
    val server = conf.getValidatedString(s"$basePath.interface")
    val port = conf.getValidatedInt(s"$basePath.port")
    val scheme = conf.getValidatedString(s"$basePath.scheme", default = Some("http"))
    (server |@| port |@| scheme) map ServiceConfig.apply
  }
}

final case class SwaggerOauthConfig(clientId: String, realm: String, appName: String)

object SwaggerOauthConfig {
  private[config] def getFromConfig(conf: Config, basePath: String): ErrorOr[SwaggerOauthConfig] = {
    def getValidatedOption(option: String) = conf.getValidatedString(s"$basePath.$option")

    (getValidatedOption("client_id") |@| getValidatedOption("realm") |@| getValidatedOption("app_name")) map SwaggerOauthConfig.apply
  }
}
