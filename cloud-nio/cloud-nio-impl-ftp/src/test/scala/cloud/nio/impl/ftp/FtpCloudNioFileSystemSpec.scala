package cloud.nio.impl.ftp

import com.typesafe.config.ConfigFactory
import io.github.andrebeat.pool.Pool.ClosedPoolException
import org.apache.commons.net.ftp.FTPClient
import org.scalatest.concurrent.Eventually
import org.scalatest.flatspec.AnyFlatSpec
import org.scalatest.matchers.should.Matchers

import scala.concurrent.TimeoutException
import scala.concurrent.duration._


class FtpCloudNioFileSystemSpec extends AnyFlatSpec with Matchers with Eventually {

  behavior of "FtpCloudNioFileSystemSpec"

  override val patienceConfig = PatienceConfig(timeout = scaled(5.seconds), interval = scaled(1.second))
  implicit val patience = patienceConfig

  it should "lease the number of clients configured, not more, not less" in {
    val fileSystems = new FtpFileSystems(FtpFileSystems.DefaultConfig.copy(leaseTimeout = Option(1.second), capacity = 3))
    val provider = new FtpCloudNioFileSystemProvider(ConfigFactory.empty, FtpAnonymousCredentials, fileSystems)
    val fileSystem = new FtpCloudNioFileSystem(provider, "ftp.example.com") {
      // Override so we don't try to connect to anything
      override private[ftp] lazy val clientFactory = () => { new FTPClient() }
    }

    val client1 = fileSystem.leaseClient
    fileSystem.leaseClient
    fileSystem.leaseClient
    // The 4th one should timeout
    a[TimeoutException] shouldBe thrownBy(fileSystem.leaseClient)
    // Release the first one and make sure it works now
    client1.release()
    fileSystem.leaseClient

    // When the filesystem is closed we can't get clients anymore
    fileSystem.close()
    a[ClosedPoolException] shouldBe thrownBy(fileSystem.leaseClient)
  }

  it should "lease a pair of clients atomically" in {
    val fileSystems = new FtpFileSystems(FtpFileSystems.DefaultConfig.copy(capacity = 2))
    val provider = new FtpCloudNioFileSystemProvider(ConfigFactory.empty, FtpAnonymousCredentials, fileSystems)
    val fileSystem = new FtpCloudNioFileSystem(provider, "ftp.example.com") {
      // Override so we don't try to connect to anything
      override private[ftp] lazy val clientFactory = () => { new FTPClient() }

      override def leaseClient = {
        val lease = super.leaseClient
        // Increases the chance that another thread will try to acquire the lease
        Thread.`yield`()
        lease
      }
    }

    val threads = (1 to 10) map { _ =>
      new Thread {
        override def run() = {
          val clients = fileSystem.leaseClientsPair
          clients._1.release()
          clients._2.release()
        }
      }
    }

    threads.foreach(_.start())

    eventually {
      threads.forall(_.getState == Thread.State.TERMINATED) shouldBe true
    }

    fileSystem.close()
  }
}
