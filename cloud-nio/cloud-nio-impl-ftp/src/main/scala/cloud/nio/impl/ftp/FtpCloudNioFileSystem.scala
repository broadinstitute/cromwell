package cloud.nio.impl.ftp

import cloud.nio.spi.CloudNioFileSystem
import com.typesafe.scalalogging.StrictLogging
import net.ceedubs.ficus.Ficus._
import org.apache.commons.net.ftp.FTPClient
import org.slf4j.LoggerFactory

import scala.concurrent.TimeoutException
import scala.concurrent.duration._

object FtpCloudNioFileSystem {
  val logger = LoggerFactory.getLogger("FtpFileSystem")
}

class FtpCloudNioFileSystem(provider: FtpCloudNioFileSystemProvider, host: String) extends CloudNioFileSystem(provider, host) with StrictLogging {
  private val credentials = provider.credentials
  private val leaseTimeout = provider.config.getAs[FiniteDuration]("acquire-connection-timeout")
  // Cannot be less than 2, otherwise we can't copy files as we need 2 connections to copy a file (one for downstream and one for upstream)
  private val capacity = Math.max(provider.config.getAs[Int]("connection-count-per-user").getOrElse(2), 2)
  private val idleConnectionTimeout = provider.config.getAs[FiniteDuration]("idle-connection-timeout").getOrElse(10.minutes)
  private val connectionPort = provider.config.getAs[Int]("connection-port").getOrElse(21)
  private val connectionModeFunction: FTPClient => Unit = provider.config.getAs[String]("connection-mode").getOrElse("passive") match {
    case "passive" => client: FTPClient => client.enterLocalPassiveMode() 
    case "active" => client: FTPClient => client.enterLocalActiveMode()
    case other => client: FTPClient => {
      logger.warn(s"Unrecognized connection mode $other, defaulting to passive mode")
      client.enterLocalPassiveMode()
    }
  }
  
  private [ftp] lazy val clientFactory = () => {
    val client = new FTPClient()
    client.setDefaultPort(connectionPort)
    client.connect(host)
    // Calling connect RESETS the mode, so make sure to call enterLocalPassiveMode AFTER connecting !!
    connectionModeFunction(client)
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

