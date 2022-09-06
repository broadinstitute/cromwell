package cromwell.filesystems.blob

import com.azure.core.credential.AzureSasCredential
import com.typesafe.config.ConfigFactory
import common.mock.MockSugar
import org.mockito.Mockito._
import org.scalatest.flatspec.AnyFlatSpec
import org.scalatest.matchers.should.Matchers

import java.nio.file.{FileSystem, FileSystemNotFoundException}
import java.time.format.DateTimeFormatter
import java.time.temporal.ChronoUnit
import java.time.{Duration, Instant, ZoneId}
import scala.util.{Failure, Try}


object BlobPathBuilderFactorySpec {
  def buildExampleSasToken(expiry: Instant): AzureSasCredential = {
    val formatter = DateTimeFormatter.ofPattern("yyyy-MM-dd").withZone(ZoneId.systemDefault())
    val sv = formatter.format(expiry)
    val se = expiry.toString().replace(":","%3A")
    new AzureSasCredential(s"sv=$sv&se=$se&sr=c&sp=rcl")
  }
}
class BlobPathBuilderFactorySpec extends AnyFlatSpec with Matchers with MockSugar {
  def generateTokenExpiration(minutes: Long) = Instant.now.plus(minutes, ChronoUnit.MINUTES)
  it should "parse configs for a functioning factory" in {
    val endpoint = BlobPathBuilderSpec.buildEndpoint("storageAccount")
    val container = BlobContainerName("storageContainer")
    val workspaceId = WorkspaceId("mockWorkspaceId")
    val workspaceManagerURL = WorkspaceManagerURL("https://test.ws.org")
    val instanceConfig = ConfigFactory.parseString(
      s"""
      |container = "$container"

      |endpoint = "$endpoint"
      |expiry-buffer-minutes = "10"
      |workspace-id = "$workspaceId"
      """.stripMargin)
    val singletonConfig = ConfigFactory.parseString(s"""workspace-manager-url = "$workspaceManagerURL" """)
    val globalConfig = ConfigFactory.parseString("""""")
    val factory = BlobPathBuilderFactory(globalConfig, instanceConfig, new BlobFileSystemConfig(singletonConfig))
    factory.container should equal(container)
    factory.endpoint should equal(endpoint)
    factory.expiryBufferMinutes should equal(10L)
    factory.workspaceId should contain(workspaceId)
    factory.workspaceManagerURL should contain(workspaceManagerURL)
  }

  it should "build an example sas token of the correct format" in {
    val testToken = BlobPathBuilderFactorySpec.buildExampleSasToken(Instant.ofEpochMilli(1603794041000L))
    val sourceToken = "sv=2020-10-27&se=2020-10-27T10%3A20%3A41Z&sr=c&sp=rcl"
    testToken.getSignature should equal(sourceToken)
  }

  it should "parse an expiration time from a sas token" in {
    val expiryTime = generateTokenExpiration(20L)
    val sasToken = BlobPathBuilderFactorySpec.buildExampleSasToken(expiryTime)
    val expiry = BlobFileSystemManager.parseTokenExpiry(sasToken)
    expiry should contain(expiryTime)
  }

  it should "verify an unexpired token will be processed as unexpired" in {
    val expiryTime = generateTokenExpiration(11L)
    val expired = BlobFileSystemManager.hasTokenExpired(expiryTime, Duration.ofMinutes(10L))
    expired shouldBe false
  }

  it should "test an expired token will be processed as expired" in {
    val expiryTime = generateTokenExpiration(9L)
    val expired = BlobFileSystemManager.hasTokenExpired(expiryTime, Duration.ofMinutes(10L))
    expired shouldBe true
  }

  it should "test that a filesystem gets closed correctly" in {
    val endpoint = BlobPathBuilderSpec.buildEndpoint("storageAccount")
    val azureUri = BlobFileSystemManager.uri(endpoint)
    val fileSystems = mock[FileSystemAPI]
    val fileSystem = mock[FileSystem]
    when(fileSystems.getFileSystem(azureUri)).thenReturn(Try(fileSystem))
    when(fileSystems.closeFileSystem(azureUri)).thenCallRealMethod()

    fileSystems.closeFileSystem(azureUri)
    verify(fileSystem, times(1)).close()
  }

  it should "test retrieveFileSystem with expired filesystem" in {
    val endpoint = BlobPathBuilderSpec.buildEndpoint("storageAccount")
    val expiredToken = generateTokenExpiration(9L)
    val refreshedToken = generateTokenExpiration(69L)
    val sasToken = BlobPathBuilderFactorySpec.buildExampleSasToken(refreshedToken)
    val container = BlobContainerName("storageContainer")
    val configMap = BlobFileSystemManager.buildConfigMap(sasToken, container)
    val azureUri = BlobFileSystemManager.uri(endpoint)

    val fileSystems = mock[FileSystemAPI]
    val blobTokenGenerator = mock[BlobTokenGenerator]
    when(blobTokenGenerator.generateAccessToken).thenReturn(Try(sasToken))

    val fsm = BlobFileSystemManager(container, endpoint, 10L, blobTokenGenerator, fileSystems, Some(expiredToken))
    fsm.getExpiry should contain(expiredToken)
    fsm.isTokenExpired shouldBe true
    fsm.retrieveFilesystem()

    fsm.getExpiry should contain(refreshedToken)
    fsm.isTokenExpired shouldBe false
    verify(fileSystems, never()).getFileSystem(azureUri)
    verify(fileSystems, times(1)).newFileSystem(azureUri, configMap)
    verify(fileSystems, times(1)).closeFileSystem(azureUri)
  }

  it should "test retrieveFileSystem with an unexpired fileSystem" in {
    val endpoint = BlobPathBuilderSpec.buildEndpoint("storageAccount")
    val initialToken = generateTokenExpiration(11L)
    val refreshedToken = generateTokenExpiration(71L)
    val sasToken = BlobPathBuilderFactorySpec.buildExampleSasToken(refreshedToken)
    val container = BlobContainerName("storageContainer")
    val configMap = BlobFileSystemManager.buildConfigMap(sasToken, container)
    val azureUri = BlobFileSystemManager.uri(endpoint)
    // Need a fake filesystem to supply the getFileSystem simulated try
    val dummyFileSystem = null

    val fileSystems = mock[FileSystemAPI]
    when(fileSystems.getFileSystem(azureUri)).thenReturn(Try(dummyFileSystem))

    val blobTokenGenerator = mock[BlobTokenGenerator]
    when(blobTokenGenerator.generateAccessToken).thenReturn(Try(sasToken))

    val fsm = BlobFileSystemManager(container, endpoint, 10L, blobTokenGenerator, fileSystems, Some(initialToken))
    fsm.getExpiry should contain(initialToken)
    fsm.isTokenExpired shouldBe false
    fsm.retrieveFilesystem()

    fsm.getExpiry should contain(initialToken)
    fsm.isTokenExpired shouldBe false
    verify(fileSystems, times(1)).getFileSystem(azureUri)
    verify(fileSystems, never()).newFileSystem(azureUri, configMap)
    verify(fileSystems, never()).closeFileSystem(azureUri)
  }

  it should "test retrieveFileSystem with an uninitialized filesystem" in {
    val endpoint = BlobPathBuilderSpec.buildEndpoint("storageAccount")
    val refreshedToken = generateTokenExpiration(71L)
    val sasToken = BlobPathBuilderFactorySpec.buildExampleSasToken(refreshedToken)
    val container = BlobContainerName("storageContainer")
    val configMap = BlobFileSystemManager.buildConfigMap(sasToken, container)
    val azureUri = BlobFileSystemManager.uri(endpoint)

    val fileSystems = mock[FileSystemAPI]
    when(fileSystems.getFileSystem(azureUri)).thenReturn(Failure(new FileSystemNotFoundException))
    val blobTokenGenerator = mock[BlobTokenGenerator]
    when(blobTokenGenerator.generateAccessToken).thenReturn(Try(sasToken))

    val fsm = BlobFileSystemManager(container, endpoint, 10L, blobTokenGenerator, fileSystems, Some(refreshedToken))
    fsm.getExpiry.isDefined shouldBe true
    fsm.isTokenExpired shouldBe false
    fsm.retrieveFilesystem()

    fsm.getExpiry should contain(refreshedToken)
    fsm.isTokenExpired shouldBe false
    verify(fileSystems, times(1)).getFileSystem(azureUri)
    verify(fileSystems, times(1)).newFileSystem(azureUri, configMap)
    verify(fileSystems, never()).closeFileSystem(azureUri)
  }

  it should "test retrieveFileSystem with an unknown filesystem" in {
    val endpoint = BlobPathBuilderSpec.buildEndpoint("storageAccount")
    val refreshedToken = generateTokenExpiration(71L)
    val sasToken = BlobPathBuilderFactorySpec.buildExampleSasToken(refreshedToken)
    val container = BlobContainerName("storageContainer")
    val configMap = BlobFileSystemManager.buildConfigMap(sasToken, container)
    val azureUri = BlobFileSystemManager.uri(endpoint)

    val fileSystems = mock[FileSystemAPI]
    val blobTokenGenerator = mock[BlobTokenGenerator]
    when(blobTokenGenerator.generateAccessToken).thenReturn(Try(sasToken))

    val fsm = BlobFileSystemManager(container, endpoint, 10L, blobTokenGenerator, fileSystems)
    fsm.getExpiry.isDefined shouldBe false
    fsm.isTokenExpired shouldBe false
    fsm.retrieveFilesystem()

    fsm.getExpiry should contain(refreshedToken)
    fsm.isTokenExpired shouldBe false
    verify(fileSystems, never()).getFileSystem(azureUri)
    verify(fileSystems, times(1)).newFileSystem(azureUri, configMap)
    verify(fileSystems, times(1)).closeFileSystem(azureUri)
  }
}
