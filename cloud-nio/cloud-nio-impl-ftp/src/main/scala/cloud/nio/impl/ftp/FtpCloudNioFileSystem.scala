package cloud.nio.impl.ftp

import cloud.nio.spi.CloudNioFileSystem
import net.ceedubs.ficus.Ficus._
import org.apache.commons.net.ftp.FTPClient
import org.slf4j.LoggerFactory

import scala.concurrent.TimeoutException
import scala.concurrent.duration._

object FtpCloudNioFileSystem {
  val logger = LoggerFactory.getLogger("FtpFileSystem")
}

class FtpCloudNioFileSystem(provider: FtpCloudNioFileSystemProvider, host: String) extends CloudNioFileSystem(provider, host) {
  private val credentials = provider.credentials
  private val leaseTimeout = provider.config.getAs[FiniteDuration]("acquire-connection-timeout")
  // Cannot be less than 2, otherwise we can't copy files as we need 2 connections to copy a file (one for downstream and one for upstream)
  private val capacity = Math.max(provider.config.getAs[Int]("connection-count-per-user").getOrElse(2), 2)
  private val idleConnectionTimeout = provider.config.getAs[FiniteDuration]("idle-connection-timeout").getOrElse(10.minutes)
  
  private [ftp] lazy val clientFactory = () => {
    val client = new FTPClient()
    client.enterLocalPassiveMode()
    client.connect(host)
    credentials.login(client)
    client
  }
  
  private val clientPool = new FtpClientPool(capacity, idleConnectionTimeout, clientFactory)
  
  def leaseClient = leaseTimeout match {
    case Some(timeout) => clientPool.tryAcquire(timeout).getOrElse(throw new TimeoutException("Timed out waiting for an available connection, try again later."))
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

