package cromiam.server.config

import akka.http.scaladsl.settings.ServerSettings
import cats.instances.option._
import cats.syntax.apply._
import cats.syntax.functor._
import cats.syntax.validated._
import com.typesafe.config.Config
import common.validation.ErrorOr.ErrorOr
import common.validation.Validation._
import cromiam.server.config.CromIamServerConfig._
import net.ceedubs.ficus.Ficus._

import scala.collection.JavaConverters._
import scala.util.{Failure, Success, Try}

final case class CromIamServerConfig(cromIamConfig: CromIamConfig,
                                     cromwellConfig: ServiceConfig,
                                     cromwellAbortConfig: ServiceConfig,
                                     samConfig: SamClientConfig,
                                     swaggerOauthConfig: SwaggerOauthConfig)

object CromIamServerConfig {
  def getFromConfig(conf: Config): ErrorOr[CromIamServerConfig] = {
    val cromIamConfig = CromIamConfig.getFromConfig(conf, "cromiam")
    val cromwellConfig = ServiceConfig.getFromConfig(conf, "cromwell")
    val cromwellAbortConfig = conf.as[Option[Config]]("cromwell_abort") as { ServiceConfig.getFromConfig(conf, "cromwell_abort") }
    val effectiveCromwellAbortConfig = cromwellAbortConfig.getOrElse(cromwellConfig)
    val samConfig = SamClientConfig.getFromConfig(conf, "sam")
    val googleConfig = SwaggerOauthConfig.getFromConfig(conf, "swagger_oauth")

    (cromIamConfig, cromwellConfig, effectiveCromwellAbortConfig, samConfig, googleConfig) mapN CromIamServerConfig.apply
  }

  private[config] def getValidatedConfigPath[A](conf: Config, path: String, getter: (Config, String) => A, default: Option[A] = None): ErrorOr[A] = {
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
    def getValidatedString(path: String, default: Option[String] = None): ErrorOr[String] = getValidatedConfigPath(conf, path, (c, p) => c.getString(p), default)
    def getValidatedInt(path: String): ErrorOr[Int] = getValidatedConfigPath(conf, path, (c, p) => c.getInt(p))
    def getValidatedStringList(path: String): ErrorOr[List[String]] = getValidatedConfigPath[List[String]](conf, path, (c, p) => c.getStringList(p).asScala.toList)
  }
}

final case class CromIamConfig(http: ServiceConfig, serverSettings: ServerSettings)

object CromIamConfig {

  private def getValidatedServerSettings(conf: Config): ErrorOr[ServerSettings] = {
    Try(ServerSettings(conf)) match {
      case Success(serverSettings) => serverSettings.validNel
      case Failure(e) =>
        s"Unable to generate server settings from configuration file: ${e.getMessage}".invalidNel
    }
  }

  private[config] def getFromConfig(conf: Config, basePath: String): ErrorOr[CromIamConfig] = {
    val serviceConfig = ServiceConfig.getFromConfig(conf, basePath)
    val serverSettings = getValidatedServerSettings(conf)

    (serviceConfig, serverSettings) mapN CromIamConfig.apply
  }
}

final case class SamClientConfig(http: ServiceConfig, checkSubmitWhitelist: Boolean)

object SamClientConfig {
  private[config] def getFromConfig(conf: Config, basePath: String): ErrorOr[SamClientConfig] = {
    val serviceConfig = ServiceConfig.getFromConfig(conf, basePath)
    val checkSubmitWhitelist = validate(conf.getOrElse(s"$basePath.check-submit-whitelist", true))

    (serviceConfig, checkSubmitWhitelist) mapN SamClientConfig.apply
  }
}

final case class ServiceConfig(interface: String, port: Int, scheme: String)

object ServiceConfig {
  private[config] def getFromConfig(conf: Config, basePath: String): ErrorOr[ServiceConfig] = {
    val server = conf.getValidatedString(s"$basePath.interface")
    val port = conf.getValidatedInt(s"$basePath.port")
    val scheme = conf.getValidatedString(s"$basePath.scheme", default = Some("http"))
    (server, port, scheme) mapN ServiceConfig.apply
  }
}

final case class SwaggerOauthConfig(clientId: String, realm: String, appName: String)

object SwaggerOauthConfig {
  private[config] def getFromConfig(conf: Config, basePath: String): ErrorOr[SwaggerOauthConfig] = {
    def getValidatedOption(option: String) = conf.getValidatedString(s"$basePath.$option")

    (getValidatedOption("client_id"), getValidatedOption("realm"), getValidatedOption("app_name")) mapN SwaggerOauthConfig.apply
   }
}
