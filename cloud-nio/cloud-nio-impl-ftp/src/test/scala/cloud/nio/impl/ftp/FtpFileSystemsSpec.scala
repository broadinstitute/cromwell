package cloud.nio.impl.ftp

import cloud.nio.impl.ftp.FtpFileSystems.FtpCacheKey
import com.typesafe.config.ConfigFactory
import common.assertion.CromwellTimeoutSpec
import org.scalamock.function.MockFunction1
import org.scalamock.scalatest.MockFactory
import org.scalatest.flatspec.AnyFlatSpec
import org.scalatest.matchers.should.Matchers

class FtpFileSystemsSpec extends AnyFlatSpec with CromwellTimeoutSpec with Matchers with MockFactory {

  behavior of "FtpFileSystemsSpec"

  val emptyConfig = ConfigFactory.empty
  
  it should "cache file systems per server per user" in {
    val mockCreateFunction = mockFunction[FtpCacheKey, FtpCloudNioFileSystem]
    val ftpFileSystems = new MockFtpFileSystems(FtpFileSystems.DefaultConfig, mockCreateFunction)
    
    val authenticatedCredentials1 = FtpAuthenticatedCredentials("user1", "password", None)
    
    val provider = new FtpCloudNioFileSystemProvider(emptyConfig, authenticatedCredentials1, ftpFileSystems)
    // Same as provider1, just other instance
    val providerClone = new FtpCloudNioFileSystemProvider(emptyConfig, authenticatedCredentials1, ftpFileSystems)
    
    // Expects the creation function to only be called once, since the 2 providers have the same credentials, even though they're
    // different instances
    mockCreateFunction.expects(FtpCacheKey("host1.com", provider))
      .returns(new FtpCloudNioFileSystem(provider, "host1.com"))
      .once()

    provider.newCloudNioFileSystemFromHost("host1.com")
    providerClone.newCloudNioFileSystemFromHost("host1.com")
  }
  
  class MockFtpFileSystems(conf: FtpFileSystemsConfiguration, mockCreateFunction: MockFunction1[FtpCacheKey, FtpCloudNioFileSystem]) extends FtpFileSystems(conf) {
    override private[ftp] def createFileSystem(key: FtpCacheKey) = mockCreateFunction(key)
  }
}
