package cloud.nio.impl.ftp

import java.net.URI
import java.nio.channels.ReadableByteChannel
import java.nio.file.FileAlreadyExistsException

import cloud.nio.impl.ftp.FtpUtil.FtpIoException
import cloud.nio.spi.{CloudNioRegularFileAttributes, CloudNioRetry}
import com.typesafe.config.ConfigFactory
import org.apache.commons.net.ftp.FTPReply
import org.scalamock.scalatest.{MixedMockFactory, MockFactory}
import org.scalatest.{FlatSpec, Matchers}

class FtpCloudNioFileSystemProviderSpec extends FlatSpec with Matchers with MockFactory with MixedMockFactory with MockFtpFileSystem {

  behavior of "FtpCloudNioFileSystemProviderSpec"

  it should "return true for isFatal if the underlying FTP exception is non retryable" in {
    mockProvider.isFatal(FtpIoException("something fatal", FTPReply.FAILED_SECURITY_CHECK, "blah")) shouldBe true
  }

  it should "return false for isFatal if the underlying FTP exception is retryable" in {
    mockProvider.isFatal(FtpIoException("something fatal", FTPReply.FILE_ACTION_NOT_TAKEN, "blah")) shouldBe false
  }

  it should "return true for isFatal if the exception is unknown" in {
    mockProvider.isFatal(new Exception) shouldBe true
  }

  it should "return false for isTransient" in {
    mockProvider.isTransient(new Exception) shouldBe false
  }

  it should "have the right scheme" in {
    mockProvider.getScheme shouldBe "ftp"
  }

  it should "not use pseudo directories" in {
    mockProvider.usePseudoDirectories shouldBe false
  }

  it should "pre compute the size before opening a read channel to avoid deadlocks" in {
    val mockSizeFunction = mockFunction[Long]
    val provider = new FtpCloudNioFileSystemProvider(ConfigFactory.empty(), FtpAnonymousCredentials, ftpFileSystems) {

      override def fileProvider = new FtpCloudNioFileProvider(this) {
        override def fileAttributes(cloudHost: String, cloudPath: String) =
          Option(
            new CloudNioRegularFileAttributes {
              override def fileHash = throw new UnsupportedOperationException()
              override def lastModifiedTime() = throw new UnsupportedOperationException()
              override def size() = mockSizeFunction()
              override def fileKey() = throw new UnsupportedOperationException()
            }
          )

        override def read(cloudHost: String, cloudPath: String, offset: Long) = {
          mock[ReadableByteChannel]
        }
      }
    }

    // This should only be called once, not every time we ask for the channel size
    mockSizeFunction.expects().onCall(_ => 60).once()
    val cloudNioPath = provider.getPath(URI.create("ftp://host.com/my_file.txt"))
    val channel = provider.cloudNioReadChannel(new CloudNioRetry(ConfigFactory.empty()), cloudNioPath)

    channel.size() shouldBe 60L
    channel.size() shouldBe 60L
  }

  it should "create a directory" in {
    val directoryPath = "/root/new_directory"
    val newDirectory = mockProvider.getPath(URI.create("ftp://localhost/root/new_directory"))

    fakeUnixFileSystem.exists(directoryPath) shouldBe false
    mockProvider.createDirectory(newDirectory)
    fakeUnixFileSystem.exists(directoryPath) shouldBe true
    
    // Now we should throw an exception because the directory exists
    a[FileAlreadyExistsException] shouldBe thrownBy(mockProvider.createDirectory(newDirectory))
  }
}
