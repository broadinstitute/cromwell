package cloud.nio.impl.ftp

import java.net.URI
import java.nio.channels.ReadableByteChannel
import java.nio.file.FileAlreadyExistsException
import cloud.nio.impl.ftp.FtpUtil.FtpIoException
import cloud.nio.spi.{CloudNioRegularFileAttributes, CloudNioRetry}
import com.typesafe.config.ConfigFactory
import common.assertion.CromwellTimeoutSpec
import common.mock.MockSugar
import org.apache.commons.net.ftp.FTPReply
import org.mockito.Mockito._
import org.scalatest.flatspec.AnyFlatSpec
import org.scalatest.matchers.should.Matchers

class FtpCloudNioFileSystemProviderSpec
    extends AnyFlatSpec
    with CromwellTimeoutSpec
    with Matchers
    with MockSugar
    with MockFtpFileSystem {

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
    val mockSizeFunction = mock[() => Long]
    val provider: FtpCloudNioFileSystemProvider = new FtpCloudNioFileSystemProvider(
      ConfigFactory.empty,
      FtpAnonymousCredentials,
      ftpFileSystems
    ) {

      override def fileProvider: FtpCloudNioFileProvider = new FtpCloudNioFileProvider(this) {
        override def fileAttributes(cloudHost: String, cloudPath: String): Option[CloudNioRegularFileAttributes] =
          Option(
            new CloudNioRegularFileAttributes {
              override def fileHashes: Map[String, String] = throw new UnsupportedOperationException()
              override def lastModifiedTime() = throw new UnsupportedOperationException()
              override def size(): Long = mockSizeFunction()
              override def fileKey() = throw new UnsupportedOperationException()
            }
          )

        override def read(cloudHost: String, cloudPath: String, offset: Long): ReadableByteChannel =
          mock[ReadableByteChannel]
      }
    }

    // This should only be called once, not every time we ask for the channel size
    when(mockSizeFunction.apply()).thenReturn(60)
    val cloudNioPath = provider.getPath(URI.create("ftp://host.com/my_file.txt"))
    val channel = provider.cloudNioReadChannel(new CloudNioRetry(ConfigFactory.empty()), cloudNioPath)
    verify(mockSizeFunction).apply()

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
