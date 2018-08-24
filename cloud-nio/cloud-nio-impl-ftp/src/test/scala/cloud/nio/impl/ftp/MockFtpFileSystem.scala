package cloud.nio.impl.ftp

import com.typesafe.config.ConfigFactory
import org.mockftpserver.fake.filesystem.{DirectoryEntry, UnixFakeFileSystem}
import org.mockftpserver.fake.{FakeFtpServer, UserAccount}
import org.scalatest.{BeforeAndAfterAll, Suite}

import scala.collection.JavaConverters._

trait MockFtpFileSystem extends BeforeAndAfterAll { this: Suite =>
  private var connectionPort: Option[Int] = None
  
  val fakeFtpServer = new FakeFtpServer()
  fakeFtpServer.setServerControlPort(0)
  fakeFtpServer.addUserAccount(new UserAccount("test_user", "test_password", "/"))
  val fakeUnixFileSystem = new UnixFakeFileSystem()
  fakeUnixFileSystem.setCreateParentDirectoriesAutomatically(false)
  fakeUnixFileSystem.add(new DirectoryEntry("/"))
  fakeUnixFileSystem.add(new DirectoryEntry("/root"))
  fakeFtpServer.setFileSystem(fakeUnixFileSystem)

  override def beforeAll = {
    fakeFtpServer.start()
    connectionPort = Option(fakeFtpServer.getServerControlPort)
  }

  override def afterAll = {
    fakeFtpServer.stop()
  }
  
  // Do not call this before starting the server
  lazy val mockProvider = {
    val config = ConfigFactory.parseMap(Map("connection-port" -> connectionPort.getOrElse(throw new RuntimeException("Fake FTP server has not been started"))).asJava)
    new FtpCloudNioFileSystemProvider(config, FtpAuthenticatedCredentials("test_user", "test_password", None))
  }
}
