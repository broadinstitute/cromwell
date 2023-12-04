package cromwell.filesystems.ftp

import cats.syntax.apply._
import cats.syntax.validated._
import cloud.nio.impl.ftp.{FtpAnonymousCredentials, FtpAuthenticatedCredentials, FtpCredentials}
import com.typesafe.config.Config
import common.validation.ErrorOr.ErrorOr
import common.validation.Validation._
import net.ceedubs.ficus.Ficus._

/**
  * Configuration for an instance of an FTPPathBuilderFactory
  * 
  * e.g:
  * 
  * engine.filesystems.ftp.config = {
  *   auth {
  *     username = "username"
  *     password = "password"
  *     account = "account"
  *   }
  * }
  */
case class FtpInstanceConfiguration(ftpCredentials: FtpCredentials)

object FtpInstanceConfiguration {
  lazy val Default = FtpInstanceConfiguration(FtpAnonymousCredentials)

  def apply(conf: Config): FtpInstanceConfiguration = {
    val credentials: ErrorOr[FtpCredentials] = conf.getAs[Config]("auth") map { authConfig =>
      val username = validate(authConfig.as[String]("username"))
      val password = validate(authConfig.as[String]("password"))
      val account = validate(authConfig.getAs[String]("account"))
      (username, password, account) mapN FtpAuthenticatedCredentials.apply
    } getOrElse Default.ftpCredentials.validNel

    val validatedConfiguration = credentials.map(FtpInstanceConfiguration.apply)

    validatedConfiguration.unsafe("FTP configuration is not valid")
  }
}
