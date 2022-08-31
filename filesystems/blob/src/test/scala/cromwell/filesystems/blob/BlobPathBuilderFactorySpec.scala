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
    val endpoint = BlobPathBuilderSpec.buildEndpoint("coaexternalstorage")
    val store = BlobContainerName("inputs")
    val sasToken = "{SAS TOKEN HERE}"
    val workspaceId = WorkspaceId("mockWorkspaceId")
    val workspaceManagerURL = WorkspaceManagerURL("https://test.ws.org")
    val instanceConfig = ConfigFactory.parseString(
      s"""
      |sas-token = "$sasToken"
      |store = "$store"
      |endpoint = "$endpoint"
      |workspace-id = "$workspaceId"
      """.stripMargin)
    val singletonConfig = ConfigFactory.parseString(s"""workspace-manager-url = "$workspaceManagerURL" """)
    val globalConfig = ConfigFactory.parseString("""""")
    val factory = BlobPathBuilderFactory(globalConfig, instanceConfig, new BlobFileSystemConfig(singletonConfig))
    factory.container should equal(store)
    factory.endpoint should equal(endpoint)
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

  it should "test retrieveFileSystem with expired filesystem" in {
    val endpoint = BlobPathBuilderSpec.buildEndpoint("coaexternalstorage")
    val expiredToken = generateTokenExpiration(9L)
    val refreshedToken = generateTokenExpiration(69L)
    val sasToken = BlobPathBuilderFactorySpec.buildExampleSasToken(refreshedToken)
    val store = BlobContainerName("inputs")
    val configMap = BlobFileSystemManager.buildConfigMap(sasToken, store)
    val azureUri = BlobFileSystemManager.uri(endpoint)

    val fileSystems = mock[FileSystemAPI]
    val fileSystem = mock[FileSystem]
    when(fileSystems.getFileSystem(azureUri)).thenReturn(Try(fileSystem))
    when(fileSystems.newFileSystem(azureUri, configMap)).thenReturn(fileSystem)
    val blobTokenGenerator = mock[BlobTokenGenerator]
    when(blobTokenGenerator.generateAccessToken).thenReturn(Try(sasToken))

    val fsm = BlobFileSystemManager(store, endpoint, 10L, blobTokenGenerator, fileSystems, Some(expiredToken))
    fsm.getExpiry() should contain(expiredToken)
    fsm.hasTokenExpired shouldBe true
    fsm.retrieveFilesystem()

    fsm.getExpiry().isDefined shouldBe true
    fsm.getExpiry() should contain(refreshedToken)
    fsm.hasTokenExpired shouldBe false
    verify(fileSystems, times(1)).getFileSystem(azureUri)
    verify(fileSystems, times(1)).newFileSystem(azureUri, configMap)
    verify(fileSystem, times(1)).close()
  }

  it should "test retrieveFileSystem with an unexpired fileSystem" in {
    val endpoint = BlobPathBuilderSpec.buildEndpoint("coaexternalstorage")
    val initialToken = generateTokenExpiration(11L)
    val refreshedToken = generateTokenExpiration(71L)
    val sasToken = BlobPathBuilderFactorySpec.buildExampleSasToken(refreshedToken)
    val store = BlobContainerName("inputs")
    val configMap = BlobFileSystemManager.buildConfigMap(sasToken, store)
    val azureUri = BlobFileSystemManager.uri(endpoint)

    val fileSystems = mock[FileSystemAPI]
    val fileSystem = mock[FileSystem]
    when(fileSystems.getFileSystem(azureUri)).thenReturn(Try(fileSystem))
    when(fileSystems.newFileSystem(azureUri, configMap)).thenReturn(fileSystem)
    val blobTokenGenerator = mock[BlobTokenGenerator]
    when(blobTokenGenerator.generateAccessToken).thenReturn(Try(sasToken))

    val fsm = BlobFileSystemManager(store, endpoint, 10L, blobTokenGenerator, fileSystems, Some(initialToken))
    fsm.getExpiry() should contain(initialToken)
    fsm.hasTokenExpired shouldBe false
    fsm.retrieveFilesystem()

    fsm.getExpiry().isDefined shouldBe true
    fsm.getExpiry() should contain(initialToken)
    fsm.hasTokenExpired shouldBe false
    verify(fileSystems, times(1)).getFileSystem(azureUri)
    verify(fileSystems, never()).newFileSystem(azureUri, configMap)
    verify(fileSystem, never()).close()
  }

  it should "test retrieveFileSystem with an uninitialized filesystem" in {
    val endpoint = BlobPathBuilderSpec.buildEndpoint("coaexternalstorage")
    val refreshedToken = generateTokenExpiration(71L)
    val sasToken = BlobPathBuilderFactorySpec.buildExampleSasToken(refreshedToken)
    val store = BlobContainerName("inputs")
    val configMap = BlobFileSystemManager.buildConfigMap(sasToken, store)
    val azureUri = BlobFileSystemManager.uri(endpoint)

    val fileSystems = mock[FileSystemAPI]
    val fileSystem = mock[FileSystem]
    when(fileSystems.getFileSystem(azureUri)).thenReturn(Failure(new FileSystemNotFoundException))
    when(fileSystems.newFileSystem(azureUri, configMap)).thenReturn(fileSystem)
    val blobTokenGenerator = mock[BlobTokenGenerator]
    when(blobTokenGenerator.generateAccessToken).thenReturn(Try(sasToken))

    val fsm = BlobFileSystemManager(store, endpoint, 10L, blobTokenGenerator, fileSystems)
    fsm.getExpiry().isDefined shouldBe false
    fsm.hasTokenExpired shouldBe false
    fsm.retrieveFilesystem()

    fsm.getExpiry().isDefined shouldBe true
    fsm.getExpiry() should contain(refreshedToken)
    fsm.hasTokenExpired shouldBe false
    verify(fileSystems, times(1)).getFileSystem(azureUri)
    verify(fileSystems, times(1)).newFileSystem(azureUri, configMap)
    verify(fileSystem, never()).close()
  }

  it should "test retrieveFileSystem with an unknown filesystem" in {
    val endpoint = BlobPathBuilderSpec.buildEndpoint("coaexternalstorage")
    val refreshedToken = generateTokenExpiration(71L)
    val sasToken = BlobPathBuilderFactorySpec.buildExampleSasToken(refreshedToken)
    val store = BlobContainerName("inputs")
    val configMap = BlobFileSystemManager.buildConfigMap(sasToken, store)
    val azureUri = BlobFileSystemManager.uri(endpoint)

    val fileSystems = mock[FileSystemAPI]
    val fileSystem = mock[FileSystem]
    when(fileSystems.getFileSystem(azureUri)).thenReturn(Try(fileSystem))
    when(fileSystems.newFileSystem(azureUri, configMap)).thenReturn(fileSystem)
    val blobTokenGenerator = mock[BlobTokenGenerator]
    when(blobTokenGenerator.generateAccessToken).thenReturn(Try(sasToken))

    val fsm = BlobFileSystemManager(store, endpoint, 10L, blobTokenGenerator, fileSystems)
    fsm.getExpiry().isDefined shouldBe false
    fsm.hasTokenExpired shouldBe false
    fsm.retrieveFilesystem()

    fsm.getExpiry().isDefined shouldBe true
    fsm.getExpiry() should contain(refreshedToken)
    fsm.hasTokenExpired shouldBe false
    verify(fileSystems, times(1)).getFileSystem(azureUri)
    verify(fileSystems, times(1)).newFileSystem(azureUri, configMap)
    verify(fileSystem, times(1)).close()
  }
}
