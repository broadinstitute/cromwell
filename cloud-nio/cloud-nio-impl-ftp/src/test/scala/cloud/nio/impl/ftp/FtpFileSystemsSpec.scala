package cloud.nio.impl.ftp

import cloud.nio.impl.ftp.FtpFileSystems.FtpCacheKey
import com.typesafe.config.ConfigFactory
import common.assertion.CromwellTimeoutSpec
import common.mock.MockSugar
import org.mockito.Mockito._
import org.scalatest.flatspec.AnyFlatSpec
import org.scalatest.matchers.should.Matchers

class FtpFileSystemsSpec extends AnyFlatSpec with CromwellTimeoutSpec with Matchers with MockSugar {

  behavior of "FtpFileSystemsSpec"

  private val emptyConfig = ConfigFactory.empty

  it should "cache file systems per server per user" in {
    val mockCreateFunction = mock[FtpCacheKey => FtpCloudNioFileSystem]
    val ftpFileSystems = new MockFtpFileSystems(FtpFileSystems.DefaultConfig, mockCreateFunction)

    val authenticatedCredentials1 = FtpAuthenticatedCredentials("user1", "password", None)

    val provider = new FtpCloudNioFileSystemProvider(emptyConfig, authenticatedCredentials1, ftpFileSystems)
    // Same as provider1, just other instance
    val providerClone = new FtpCloudNioFileSystemProvider(emptyConfig, authenticatedCredentials1, ftpFileSystems)

    // Expects the creation function to only be called once, since the 2 providers have the same credentials, even though they're
    // different instances
    val ftpCacheKey = FtpCacheKey("host1.com", provider)
    when(mockCreateFunction.apply(ftpCacheKey)) thenReturn
      new FtpCloudNioFileSystem(provider, "host1.com")

    provider.newCloudNioFileSystemFromHost("host1.com")
    providerClone.newCloudNioFileSystemFromHost("host1.com")
    verify(mockCreateFunction).apply(ftpCacheKey)
  }

  class MockFtpFileSystems(conf: FtpFileSystemsConfiguration,
                           mockCreateFunction: FtpCacheKey => FtpCloudNioFileSystem) extends FtpFileSystems(conf) {
    override private[ftp] def createFileSystem(key: FtpCacheKey) = mockCreateFunction(key)
  }
}
