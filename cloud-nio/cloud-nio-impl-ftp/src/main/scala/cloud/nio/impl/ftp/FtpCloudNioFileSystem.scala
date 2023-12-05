package cloud.nio.impl.ftp

import cloud.nio.impl.ftp.FtpFileSystemsConfiguration.{Active, Passive}
import cloud.nio.spi.CloudNioFileSystem
import com.typesafe.scalalogging.StrictLogging
import org.apache.commons.net.ftp.FTPClient
import org.slf4j.LoggerFactory

import scala.concurrent.TimeoutException

object FtpCloudNioFileSystem {
  val logger = LoggerFactory.getLogger("FtpFileSystem")
}

class FtpCloudNioFileSystem(provider: FtpCloudNioFileSystemProvider, host: String)
    extends CloudNioFileSystem(provider, host)
    with StrictLogging {
  private val credentials = provider.credentials
  private val ftpConfig = provider.ftpConfig
  private val connectionModeFunction: FTPClient => Unit = ftpConfig.connectionMode match {
    case Passive => client: FTPClient => client.enterLocalPassiveMode()
    case Active => client: FTPClient => client.enterLocalActiveMode()
  }

  private[ftp] lazy val clientFactory = () => {
    val client = new FTPClient()
    client.setDefaultPort(ftpConfig.connectionPort)
    client.connect(host)
    // Calling connect RESETS the mode, so make sure to call enterLocalPassiveMode AFTER connecting !!
    connectionModeFunction(client)
    credentials.login(client)
    client
  }

  private val clientPool = new FtpClientPool(ftpConfig.capacity, ftpConfig.idleConnectionTimeout, clientFactory)

  def leaseClient = ftpConfig.leaseTimeout match {
    case Some(timeout) =>
      clientPool
        .tryAcquire(timeout)
        .getOrElse(throw new TimeoutException("Timed out waiting for an available connection, try again later."))
    case _ => clientPool.acquire()
  }

  // Needs to be synchronized to avoid deadlocks
  def leaseClientsPair = clientPool.synchronized {
    (leaseClient, leaseClient)
  }

  override def close() = {
    clientPool.close()
    super.close()
  }
}
