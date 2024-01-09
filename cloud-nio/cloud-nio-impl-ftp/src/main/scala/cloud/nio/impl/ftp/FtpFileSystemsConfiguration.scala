package cloud.nio.impl.ftp

import cloud.nio.impl.ftp.FtpFileSystemsConfiguration.ConnectionMode

import scala.concurrent.duration.FiniteDuration

case class FtpFileSystemsConfiguration(cacheTTL: FiniteDuration,
                                       leaseTimeout: Option[FiniteDuration],
                                       capacity: Int,
                                       idleConnectionTimeout: FiniteDuration,
                                       connectionPort: Int,
                                       connectionMode: ConnectionMode
)

object FtpFileSystemsConfiguration {
  sealed trait ConnectionMode
  case object Passive extends ConnectionMode
  case object Active extends ConnectionMode
}
