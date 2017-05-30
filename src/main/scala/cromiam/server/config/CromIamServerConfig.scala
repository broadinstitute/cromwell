package cromiam.server.config

import akka.http.scaladsl.settings.ServerSettings
import cats.syntax.validated._
import cats.syntax.cartesian._
import com.typesafe.config.Config
import lenthall.validation.ErrorOr.ErrorOr

import scala.collection.JavaConverters._
import scala.util.{Failure, Success, Try}
import cromiam.server.config.CromIamServerConfig._

final case class CromIamServerConfig(cromIamConfig: CromIamConfig, cromwellConfig: ServiceConfig, samConfig: ServiceConfig)

object CromIamServerConfig {
  def getFromConfig(conf: Config): ErrorOr[CromIamServerConfig] = {
    val cromIamConfig = CromIamConfig.getFromConfig(conf, "cromiam")
    val cromwellConfig = ServiceConfig.getFromConfig(conf, "cromwell")
    val samConfig = ServiceConfig.getFromConfig(conf, "sam")

    (cromIamConfig |@| cromwellConfig |@| samConfig) map CromIamServerConfig.apply
  }

  private[config] def getValidatedConfigPath[A](typename: String, conf: Config, path: String, getter: (Config, String) => A): ErrorOr[A] = {
    if (conf.hasPath(path)) {
      Try(getter.apply(conf, path)) match {
        case Success(s) => s.validNel
        case Failure(e) => s"Unable to read valid value at '$path': ${e.getMessage}".invalidNel
      }
    } else {
      s"Configuration does not have path $path".invalidNel
    }
  }

  private[config] implicit final class ValidatingConfig(val conf: Config) extends AnyVal {
    def getValidatedString(path: String): ErrorOr[String] = getValidatedConfigPath("string", conf, path, (c, p) => c.getString(p))
    def getValidatedInt(path: String): ErrorOr[Int] = getValidatedConfigPath("integer", conf, path, (c, p) => c.getInt(p))
    def getValidatedStringList(path: String): ErrorOr[List[String]] = getValidatedConfigPath[List[String]]("string list", conf, path, (c, p) => c.getStringList(p).asScala.toList)
    def getValidatedServerSettings(path: String): ErrorOr[ServerSettings] = getValidatedConfigPath[Config]("server settings", conf, path, (c, p) => c.getConfig(p)) map { config => ServerSettings(config) }
  }
}

final case class CromIamConfig(userIdHeader: String, allowedUsers: List[String], http: ServiceConfig, serverSettings: ServerSettings)

object CromIamConfig {
  private[config] def getFromConfig(conf: Config, basePath: String): ErrorOr[CromIamConfig] = {
    val userIdHeader = conf.getValidatedString(s"$basePath.user_id_header")
    val allowedUserList = conf.getValidatedStringList(s"$basePath.allowed_users")
    val serviceConfig = ServiceConfig.getFromConfig(conf, s"$basePath.http")
    val serverSettings = conf.getValidatedServerSettings(s"$basePath.http")

    (userIdHeader |@| allowedUserList |@| serviceConfig |@| serverSettings) map CromIamConfig.apply
  }
}

final case class ServiceConfig(interface: String, port: Int)

object ServiceConfig {
  private[config] def getFromConfig(conf: Config, basePath: String): ErrorOr[ServiceConfig] = {
    val server = conf.getValidatedString(s"$basePath.interface")
    val port = conf.getValidatedInt(s"$basePath.port")
    (server |@| port) map ServiceConfig.apply
  }
}
