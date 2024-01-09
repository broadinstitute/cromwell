package cromwell.filesystems.blob

import com.azure.core.credential.AzureSasCredential
import com.azure.storage.blob.nio.AzureFileSystem
import common.mock.MockSugar
import org.mockito.Mockito._
import org.scalatest.flatspec.AnyFlatSpec
import org.scalatest.matchers.should.Matchers

import java.nio.file.FileSystemNotFoundException
import java.time.format.DateTimeFormatter
import java.time.temporal.ChronoUnit
import java.time.{Duration, Instant, ZoneId}
import java.util.UUID
import scala.util.{Failure, Success, Try}

object BlobPathBuilderFactorySpec {
  def buildExampleSasToken(expiry: Instant): AzureSasCredential = {
    val formatter = DateTimeFormatter.ofPattern("yyyy-MM-dd").withZone(ZoneId.systemDefault())
    val sv = formatter.format(expiry)
    val se = expiry.toString().replace(":", "%3A")
    new AzureSasCredential(s"sv=$sv&se=$se&sr=c&sp=rcl")
  }
}
class BlobPathBuilderFactorySpec extends AnyFlatSpec with Matchers with MockSugar {
  def generateTokenExpiration(minutes: Long) = Instant.now.plus(minutes, ChronoUnit.MINUTES)

  it should "build an example sas token of the correct format" in {
    val testToken = BlobPathBuilderFactorySpec.buildExampleSasToken(Instant.ofEpochMilli(1603794041000L))
    val sourceToken = "sv=2020-10-27&se=2020-10-27T10%3A20%3A41Z&sr=c&sp=rcl"
    testToken.getSignature should equal(sourceToken)
  }

  it should "test that a filesystem gets closed correctly" in {
    val endpoint = BlobPathBuilderSpec.buildEndpoint("storageAccount")
    val container = BlobContainerName("test")
    val azureUri = BlobFileSystemManager.combinedEnpointContainerUri(endpoint, container)
    val fileSystems = mock[AzureFileSystemAPI]
    val fileSystem = mock[AzureFileSystem]
    when(fileSystems.getFileSystem(azureUri)).thenReturn(Try(fileSystem))
    when(fileSystems.closeFileSystem(azureUri)).thenCallRealMethod()

    fileSystems.closeFileSystem(azureUri)
    verify(fileSystem, times(1)).close()
  }

  it should "test retrieveFileSystem with expired Terra filesystem" in {
    val endpoint = BlobPathBuilderSpec.buildEndpoint("storageAccount")
    // val expiredToken = generateTokenExpiration(9L)
    val refreshedToken = generateTokenExpiration(69L)
    val sasToken = BlobPathBuilderFactorySpec.buildExampleSasToken(refreshedToken)
    val container = BlobContainerName("sc-" + UUID.randomUUID().toString())
    val configMap = BlobFileSystemManager.buildConfigMap(sasToken, container)
    val azureUri = BlobFileSystemManager.combinedEnpointContainerUri(endpoint, container)

    // Mocking this final class requires the plugin Mock Maker Inline plugin, configured here
    // at filesystems/blob/src/test/resources/mockito-extensions/org.mockito.plugins.MockMaker
    val azureFileSystem = mock[AzureFileSystem]
    when(azureFileSystem.isExpired(Duration.ofMinutes(10L))).thenReturn(true)
    val fileSystems = mock[AzureFileSystemAPI]
    when(fileSystems.getFileSystem(azureUri)).thenReturn(Success(azureFileSystem))
    val blobTokenGenerator = mock[BlobSasTokenGenerator]
    when(blobTokenGenerator.generateBlobSasToken(endpoint, container)).thenReturn(Try(sasToken))

    val fsm = new BlobFileSystemManager(10L, blobTokenGenerator, fileSystems)
    fsm.retrieveFilesystem(endpoint, container)

    verify(fileSystems, times(1)).getFileSystem(azureUri)
    verify(fileSystems, times(1)).newFileSystem(azureUri, configMap)
    verify(fileSystems, times(1)).closeFileSystem(azureUri)
  }

  it should "test retrieveFileSystem with an unexpired Terra fileSystem" in {
    val endpoint = BlobPathBuilderSpec.buildEndpoint("storageAccount")
    // val initialToken = generateTokenExpiration(11L)
    val refreshedToken = generateTokenExpiration(71L)
    val sasToken = BlobPathBuilderFactorySpec.buildExampleSasToken(refreshedToken)
    val container = BlobContainerName("sc-" + UUID.randomUUID().toString())
    val configMap = BlobFileSystemManager.buildConfigMap(sasToken, container)
    val azureUri = BlobFileSystemManager.combinedEnpointContainerUri(endpoint, container)

    // Mocking this final class requires the plugin Mock Maker Inline plugin, configured here
    // at filesystems/blob/src/test/resources/mockito-extensions/org.mockito.plugins.MockMaker
    val azureFileSystem = mock[AzureFileSystem]
    when(azureFileSystem.isExpired(Duration.ofMinutes(10L))).thenReturn(false)
    val fileSystems = mock[AzureFileSystemAPI]
    when(fileSystems.getFileSystem(azureUri)).thenReturn(Try(azureFileSystem))

    val blobTokenGenerator = mock[BlobSasTokenGenerator]
    when(blobTokenGenerator.generateBlobSasToken(endpoint, container)).thenReturn(Try(sasToken))

    val fsm = new BlobFileSystemManager(10L, blobTokenGenerator, fileSystems)
    fsm.retrieveFilesystem(endpoint, container)

    verify(fileSystems, times(1)).getFileSystem(azureUri)
    verify(fileSystems, never()).newFileSystem(azureUri, configMap)
    verify(fileSystems, never()).closeFileSystem(azureUri)
  }

  it should "test retrieveFileSystem with an uninitialized Terra filesystem" in {
    val endpoint = BlobPathBuilderSpec.buildEndpoint("storageAccount")
    val refreshedToken = generateTokenExpiration(71L)
    val sasToken = BlobPathBuilderFactorySpec.buildExampleSasToken(refreshedToken)
    val container = BlobContainerName("sc-" + UUID.randomUUID().toString())
    val configMap = BlobFileSystemManager.buildConfigMap(sasToken, container)
    val azureUri = BlobFileSystemManager.combinedEnpointContainerUri(endpoint, container)

    // Mocking this final class requires the plugin Mock Maker Inline plugin, configured here
    // at filesystems/blob/src/test/resources/mockito-extensions/org.mockito.plugins.MockMaker
    val azureFileSystem = mock[AzureFileSystem]
    when(azureFileSystem.isExpired(Duration.ofMinutes(10L))).thenReturn(false)
    val fileSystems = mock[AzureFileSystemAPI]
    when(fileSystems.getFileSystem(azureUri)).thenReturn(Failure(new FileSystemNotFoundException))
    when(fileSystems.newFileSystem(azureUri, configMap)).thenReturn(Try(azureFileSystem))
    val blobTokenGenerator = mock[BlobSasTokenGenerator]
    when(blobTokenGenerator.generateBlobSasToken(endpoint, container)).thenReturn(Try(sasToken))

    val fsm = new BlobFileSystemManager(0L, blobTokenGenerator, fileSystems)
    fsm.retrieveFilesystem(endpoint, container)

    verify(fileSystems, times(1)).getFileSystem(azureUri)
    verify(fileSystems, times(1)).newFileSystem(azureUri, configMap)
    verify(fileSystems, never()).closeFileSystem(azureUri)
  }

  it should "test retrieveFileSystem with expired non-Terra filesystem" in {
    val endpoint = BlobPathBuilderSpec.buildEndpoint("storageAccount")
    val sasToken = BlobFileSystemManager.PLACEHOLDER_TOKEN
    val container = BlobContainerName("sc-" + UUID.randomUUID().toString())
    val configMap = BlobFileSystemManager.buildConfigMap(sasToken, container)
    val azureUri = BlobFileSystemManager.combinedEnpointContainerUri(endpoint, container)

    // Mocking this final class requires the plugin Mock Maker Inline plugin, configured here
    // at filesystems/blob/src/test/resources/mockito-extensions/org.mockito.plugins.MockMaker
    val azureFileSystem = mock[AzureFileSystem]
    when(azureFileSystem.isExpired(Duration.ofMinutes(10L))).thenReturn(true)
    val fileSystems = mock[AzureFileSystemAPI]
    when(fileSystems.getFileSystem(azureUri)).thenReturn(Success(azureFileSystem))
    val blobTokenGenerator = mock[BlobSasTokenGenerator]
    when(blobTokenGenerator.generateBlobSasToken(endpoint, container)).thenReturn(Try(sasToken))

    val fsm = new BlobFileSystemManager(10L, blobTokenGenerator, fileSystems)
    fsm.retrieveFilesystem(endpoint, container)

    verify(fileSystems, times(1)).getFileSystem(azureUri)
    verify(fileSystems, times(1)).newFileSystem(azureUri, configMap)
    verify(fileSystems, times(1)).closeFileSystem(azureUri)
  }

  it should "test retrieveFileSystem with an unexpired non-Terra fileSystem" in {
    val endpoint = BlobPathBuilderSpec.buildEndpoint("storageAccount")
    val sasToken = BlobFileSystemManager.PLACEHOLDER_TOKEN
    val container = BlobContainerName("sc-" + UUID.randomUUID().toString())
    val configMap = BlobFileSystemManager.buildConfigMap(sasToken, container)
    val azureUri = BlobFileSystemManager.combinedEnpointContainerUri(endpoint, container)

    // Mocking this final class requires the plugin Mock Maker Inline plugin, configured here
    // at filesystems/blob/src/test/resources/mockito-extensions/org.mockito.plugins.MockMaker
    val azureFileSystem = mock[AzureFileSystem]
    when(azureFileSystem.isExpired(Duration.ofMinutes(10L))).thenReturn(false)
    val fileSystems = mock[AzureFileSystemAPI]
    when(fileSystems.getFileSystem(azureUri)).thenReturn(Try(azureFileSystem))

    val blobTokenGenerator = mock[BlobSasTokenGenerator]
    when(blobTokenGenerator.generateBlobSasToken(endpoint, container)).thenReturn(Try(sasToken))

    val fsm = new BlobFileSystemManager(10L, blobTokenGenerator, fileSystems)
    fsm.retrieveFilesystem(endpoint, container)

    verify(fileSystems, times(1)).getFileSystem(azureUri)
    verify(fileSystems, never()).newFileSystem(azureUri, configMap)
    verify(fileSystems, never()).closeFileSystem(azureUri)
  }

  it should "test retrieveFileSystem with an uninitialized non-Terra filesystem" in {
    val endpoint = BlobPathBuilderSpec.buildEndpoint("storageAccount")
    val sasToken = BlobFileSystemManager.PLACEHOLDER_TOKEN
    val container = BlobContainerName("sc-" + UUID.randomUUID().toString())
    val configMap = BlobFileSystemManager.buildConfigMap(sasToken, container)
    val azureUri = BlobFileSystemManager.combinedEnpointContainerUri(endpoint, container)

    // Mocking this final class requires the plugin Mock Maker Inline plugin, configured here
    // at filesystems/blob/src/test/resources/mockito-extensions/org.mockito.plugins.MockMaker
    val azureFileSystem = mock[AzureFileSystem]
    when(azureFileSystem.isExpired(Duration.ofMinutes(10L))).thenReturn(false)
    val fileSystems = mock[AzureFileSystemAPI]
    when(fileSystems.getFileSystem(azureUri)).thenReturn(Failure(new FileSystemNotFoundException))
    when(fileSystems.newFileSystem(azureUri, configMap)).thenReturn(Try(azureFileSystem))
    val blobTokenGenerator = mock[BlobSasTokenGenerator]
    when(blobTokenGenerator.generateBlobSasToken(endpoint, container)).thenReturn(Try(sasToken))

    val fsm = new BlobFileSystemManager(0L, blobTokenGenerator, fileSystems)
    fsm.retrieveFilesystem(endpoint, container)

    verify(fileSystems, times(1)).getFileSystem(azureUri)
    verify(fileSystems, times(1)).newFileSystem(azureUri, configMap)
    verify(fileSystems, never()).closeFileSystem(azureUri)
  }
}
