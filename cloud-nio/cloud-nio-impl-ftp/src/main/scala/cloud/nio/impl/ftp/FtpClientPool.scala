package cloud.nio.impl.ftp

import com.typesafe.scalalogging.StrictLogging
import io.github.andrebeat.pool.{ExpiringPool, ReferenceType}
import org.apache.commons.net.ftp.FTPClient

import scala.concurrent.duration._

object FtpClientPool extends StrictLogging {
  def dispose(ftpClient: FTPClient) = try
    if (ftpClient.isConnected) {
      ftpClient.logout()
      ftpClient.disconnect()
    }
  catch {
    case e: Exception => logger.debug("Failed to disconnect ftp client", e)
  }
}

class FtpClientPool(capacity: Int, maxIdleTime: FiniteDuration, factory: () => FTPClient)
    extends ExpiringPool[FTPClient](
      capacity = capacity,
      maxIdleTime = maxIdleTime,
      referenceType = ReferenceType.Strong,
      _factory = factory,
      // Reset is called every time a client is added or released back to the pool. We don't want to actually reset the connection here
      // otherwise we'd need to login again and reconfigure the connection every time
      _reset = Function.const(()),
      _dispose = FtpClientPool.dispose,
      // Could not find a good health check at the moment (isAvailable and isConnected on the socket seem to both return false sometimes even if the client is fine)
      _healthCheck = Function.const(true)
    )
