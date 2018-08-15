package cloud.nio.impl.ftp

import cloud.nio.impl.ftp.FtpUtil.FtpOperation
import io.github.andrebeat.pool.{ExpiringPool, ReferenceType}
import org.apache.commons.net.ftp.FTPClient
import org.slf4j.LoggerFactory

import scala.concurrent.duration._

object FtpClientPool {
  case class DataConnectionToken(operation: FtpOperation)
  val logger = LoggerFactory.getLogger("FtpClient")
  
  def dispose(ftpClient: FTPClient) = try {
    if (ftpClient.isConnected) {
      ftpClient.logout()
      ftpClient.disconnect()
    }
  } catch {
    case e: Exception => logger.error("Failed to disconnect ftp client", e)
  }
}

class FtpClientPool(capacity: Int, maxIdleTime: FiniteDuration, factory: () => FTPClient) extends ExpiringPool[FTPClient](
  capacity = capacity,
  maxIdleTime = maxIdleTime,
  referenceType = ReferenceType.Strong,
  _factory = factory,
  _reset = Function.const(()),
  _dispose = FtpClientPool.dispose,
  _healthCheck = Function.const(true)
)
