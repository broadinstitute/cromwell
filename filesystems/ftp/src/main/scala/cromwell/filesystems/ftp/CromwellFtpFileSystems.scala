package cromwell.filesystems.ftp

import cats.syntax.apply._
import cloud.nio.impl.ftp.FtpFileSystemsConfiguration.{Active, ConnectionMode, Passive}
import cloud.nio.impl.ftp.{FtpFileSystems, FtpFileSystemsConfiguration}
import com.typesafe.config.Config
import com.typesafe.config.ConfigException.{BadValue, Missing}
import common.validation.ErrorOr.ErrorOr
import common.validation.Validation._
import cromwell.filesystems.ftp.CromwellFtpFileSystems._
import net.ceedubs.ficus.Ficus._
import net.ceedubs.ficus.readers.ValueReader

import scala.concurrent.duration.FiniteDuration

object CromwellFtpFileSystems {
  implicit val connectionModeReader = new ValueReader[ConnectionMode] {
    override def read(config: Config, path: String) = if (config.hasPath(path)) {
      config.as[String](path) match {
        case "passive" => Passive
        case "active" => Active
        case other => throw new BadValue(path, s"Unrecognized connection mode: $other")
      }
    } else throw new Missing(path)
  }

  def parseConfig(config: Config): FtpFileSystemsConfiguration = {
    val cacheTTL = validate[FiniteDuration](config.as[FiniteDuration]("cache-ttl"))
    val leaseTimeout = validate[Option[FiniteDuration]](config.getAs[FiniteDuration]("obtain-connection-timeout"))
    // Cannot be less than 2, otherwise we can't copy files as we need 2 connections to copy a file (one for downstream and one for upstream)
    val capacity: ErrorOr[Int] = validate[Int](config.as[Int]("max-connection-per-server-per-user")) map { c =>
      Math.max(2, c)
    }
    val idleConnectionTimeout = validate[FiniteDuration](config.as[FiniteDuration]("idle-connection-timeout"))
    val connectionPort = validate[Int](config.as[Int]("connection-port"))
    val connectionMode = validate[ConnectionMode](config.as[ConnectionMode]("connection-mode"))

    (cacheTTL, leaseTimeout, capacity, idleConnectionTimeout, connectionPort, connectionMode)
      .mapN(FtpFileSystemsConfiguration.apply)
      .unsafe("Failed to parse FTP configuration")
  }

  val Default = new CromwellFtpFileSystems(FtpFileSystems.Default)
}

/**
  * Cromwell Wrapper around an FtpFileSystems instance to handle parsing and validation of the configuration.
  * This is the class used as the global configuration class in the ftp filesystem configuration
  */
class CromwellFtpFileSystems(val ftpFileSystems: FtpFileSystems) {
  def this(config: Config) = this(new FtpFileSystems(parseConfig(config)))
}
