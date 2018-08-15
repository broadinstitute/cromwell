package cromwell.filesystems.ftp

import cats.syntax.apply._
import cats.syntax.validated._
import cloud.nio.impl.ftp.{FtpAnonymousCredentials, FtpAuthenticatedCredentials, FtpCredentials}
import com.typesafe.config.Config
import common.validation.ErrorOr.ErrorOr
import common.validation.Validation._
import net.ceedubs.ficus.Ficus._

import scala.concurrent.duration._

case class FtpConfiguration(ftpCredentials: FtpCredentials,
                            cacheTTL: FiniteDuration) {

}

object FtpConfiguration {
  lazy val Default = FtpConfiguration(FtpAnonymousCredentials, 24.hours)
  
  def apply(conf: Config): FtpConfiguration = {
    val credentials: ErrorOr[FtpCredentials] = conf.getAs[Config]("auth") map { authConfig =>
      val username = validate { authConfig.as[String]("username") }
      val password = validate { authConfig.as[String]("password") }
      val account = validate { authConfig.getAs[String]("account") }
      (username, password, account) mapN FtpAuthenticatedCredentials.apply
    } getOrElse Default.ftpCredentials.validNel

    val cacheTTL: ErrorOr[FiniteDuration] = validate { conf.getAs[FiniteDuration]("cache-ttl").getOrElse(Default.cacheTTL) }

    val validatedConfiguration = (credentials, cacheTTL) mapN FtpConfiguration.apply

    validatedConfiguration.unsafe("FTP configuration is not valid")
  }
}
