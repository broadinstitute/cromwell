package cromwell.filesystems.blob

import com.azure.core.credential.AzureSasCredential
import common.mock.MockSugar
import org.mockito.Mockito._
import org.scalatest.flatspec.AnyFlatSpec
import org.scalatest.matchers.should.Matchers

import java.nio.file.{FileSystem, FileSystemNotFoundException}
import java.time.format.DateTimeFormatter
import java.time.temporal.ChronoUnit
import java.time.{Duration, Instant, ZoneId}
import scala.util.{Failure, Success, Try}


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
    val buffer = Duration.ofMinutes(10L)

    val fileSystems = mock[FileSystemAPI]
    val blobTokenGenerator = mock[BlobSasTokenGenerator]
    when(blobTokenGenerator.findBlobSasToken(endpoint, container, buffer)).thenReturn(Try(sasToken))

    val fsm = new BlobFileSystemManager(10L, blobTokenGenerator, fileSystems, Some(expiredToken))
    fsm.getExpiry should contain(expiredToken)
    fsm.isTokenExpired shouldBe true
    fsm.retrieveFilesystem(endpoint, container)

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
    val buffer = Duration.ofMinutes(10L)
    // Need a fake filesystem to supply the getFileSystem simulated try
    val dummyFileSystem = mock[FileSystem]

    val fileSystems = mock[FileSystemAPI]
    when(fileSystems.getFileSystem(azureUri)).thenReturn(Try(dummyFileSystem))

    val blobTokenGenerator = mock[BlobSasTokenGenerator]
    when(blobTokenGenerator.findBlobSasToken(endpoint, container, buffer)).thenReturn(Try(sasToken))

    val fsm = new BlobFileSystemManager(10L, blobTokenGenerator, fileSystems, Some(initialToken), Some((endpoint, container)))
    fsm.getExpiry should contain(initialToken)
    fsm.isTokenExpired shouldBe false
    fsm.retrieveFilesystem(endpoint, container)

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
    val buffer = Duration.ofMinutes(10L)

    val fileSystems = mock[FileSystemAPI]
    when(fileSystems.getFileSystem(azureUri)).thenReturn(Failure(new FileSystemNotFoundException))
    val blobTokenGenerator = mock[BlobSasTokenGenerator]
    when(blobTokenGenerator.findBlobSasToken(endpoint, container, buffer)).thenReturn(Try(sasToken))

    val fsm = new BlobFileSystemManager(10L, blobTokenGenerator, fileSystems, Some(refreshedToken), Some((endpoint, container)))
    fsm.getExpiry.isDefined shouldBe true
    fsm.isTokenExpired shouldBe false
    fsm.retrieveFilesystem(endpoint, container)

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
    val buffer = Duration.ofMinutes(10L)

    val fileSystems = mock[FileSystemAPI]
    val blobTokenGenerator = mock[BlobSasTokenGenerator]
    when(blobTokenGenerator.findBlobSasToken(endpoint, container, buffer)).thenReturn(Try(sasToken))

    val fsm = new BlobFileSystemManager(10L, blobTokenGenerator, fileSystems)
    fsm.getExpiry.isDefined shouldBe false
    fsm.isTokenExpired shouldBe false
    fsm.retrieveFilesystem(endpoint, container)

    fsm.getExpiry should contain(refreshedToken)
    fsm.isTokenExpired shouldBe false
    verify(fileSystems, never()).getFileSystem(azureUri)
    verify(fileSystems, times(1)).newFileSystem(azureUri, configMap)
    verify(fileSystems, times(1)).closeFileSystem(azureUri)
  }

  it should "test retrieveFileSystem with a different fileSystem than the last retrieved" in {
    val endpoint = BlobPathBuilderSpec.buildEndpoint("storageAccount")
    val initialToken = generateTokenExpiration(11L)
    val refreshedToken = generateTokenExpiration(71L)
    val sasToken = BlobPathBuilderFactorySpec.buildExampleSasToken(refreshedToken)
    val container = BlobContainerName("storageContainer")
    val configMap = BlobFileSystemManager.buildConfigMap(sasToken, container)
    val azureUri = BlobFileSystemManager.uri(endpoint)
    val buffer = Duration.ofMinutes(10L)
    // Need a fake filesystem to supply the getFileSystem simulated try
    val dummyFileSystem = mock[FileSystem]

    val fileSystems = mock[FileSystemAPI]
    when(fileSystems.getFileSystem(azureUri)).thenReturn(Try(dummyFileSystem))

    val blobTokenGenerator = mock[BlobSasTokenGenerator]
    when(blobTokenGenerator.findBlobSasToken(endpoint, container, buffer)).thenReturn(Try(sasToken))

    val fsm = new BlobFileSystemManager(10L, blobTokenGenerator, fileSystems, Some(initialToken), Some((endpoint, BlobContainerName("oldStorageContainer"))))
    fsm.getExpiry should contain(initialToken)
    fsm.isTokenExpired shouldBe false
    fsm.retrieveFilesystem(endpoint, container)

    fsm.getExpiry should contain(refreshedToken)
    fsm.isTokenExpired shouldBe false
    verify(fileSystems, never()).getFileSystem(azureUri)
    verify(fileSystems, times(1)).newFileSystem(azureUri, configMap)
    verify(fileSystems, times(1)).closeFileSystem(azureUri)
  }

  it should "test use cached SAS token rather than requesting a new one" in {
    val endpoint = BlobPathBuilderSpec.buildEndpoint("storageAccount")
    val container = BlobContainerName("storageContainer")
    val initialToken = generateTokenExpiration(12L)
    val refreshedToken = generateTokenExpiration(72L)
    val buffer = Duration.ofMinutes(10L)
    val sasTokenOld = BlobPathBuilderFactorySpec.buildExampleSasToken(initialToken)
    val sasTokenNew = BlobPathBuilderFactorySpec.buildExampleSasToken(refreshedToken)
    val blobTokenGenerator = mock[WSMBlobSasTokenGenerator]
    when(blobTokenGenerator.findBlobSasToken(endpoint, container, buffer)).thenCallRealMethod()
    when(blobTokenGenerator.generateBlobSasToken(endpoint, container)).thenReturn(Try(sasTokenNew))
    when(blobTokenGenerator.getAvailableCachedSasToken(endpoint, container)).thenReturn(Some(sasTokenOld))
    when(blobTokenGenerator.putAvailableCachedSasToken(endpoint, container, sasTokenNew)).thenCallRealMethod()
    val sas: Try[AzureSasCredential] = blobTokenGenerator.findBlobSasToken(endpoint, container, buffer)
    BlobFileSystemManager.isSasValid(sasTokenOld, buffer) shouldBe(true)
    BlobFileSystemManager.isSasValid(sasTokenNew, buffer) shouldBe(true)
    verify(blobTokenGenerator, never()).generateBlobSasToken(endpoint, container)
    verify(blobTokenGenerator, never()).putAvailableCachedSasToken(endpoint, container, sasTokenNew)
    verify(blobTokenGenerator, times(1)).getAvailableCachedSasToken(endpoint, container)
    sas shouldBe Success(sasTokenOld)
  }

  it should "test requesting SAS token when cached one has expired" in {
    val endpoint = BlobPathBuilderSpec.buildEndpoint("storageAccount")
    val container = BlobContainerName("storageContainer")
    val initialToken = generateTokenExpiration(9L)
    val refreshedToken = generateTokenExpiration(69L)
    val buffer = Duration.ofMinutes(10L)
    val sasTokenOld = BlobPathBuilderFactorySpec.buildExampleSasToken(initialToken)
    val sasTokenNew = BlobPathBuilderFactorySpec.buildExampleSasToken(refreshedToken)
    val blobTokenGenerator = mock[WSMBlobSasTokenGenerator]
    when(blobTokenGenerator.findBlobSasToken(endpoint, container, buffer)).thenCallRealMethod()
    when(blobTokenGenerator.generateBlobSasToken(endpoint, container)).thenReturn(Success(sasTokenNew))
    when(blobTokenGenerator.getAvailableCachedSasToken(endpoint, container)).thenReturn(Some(sasTokenOld))
    val sas: Try[AzureSasCredential] = blobTokenGenerator.findBlobSasToken(endpoint, container, buffer)
    BlobFileSystemManager.isSasValid(sasTokenOld, buffer) shouldBe(false)
    BlobFileSystemManager.isSasValid(sasTokenNew, buffer) shouldBe(true)
    verify(blobTokenGenerator, times(1)).generateBlobSasToken(endpoint, container)
    verify(blobTokenGenerator, times(1)).putAvailableCachedSasToken(endpoint, container, sasTokenNew)
    verify(blobTokenGenerator, times(1)).getAvailableCachedSasToken(endpoint, container)
    sas shouldBe Success(sasTokenNew)
  }

  it should "test requesting SAS token when no cached value is available" in {
    val endpoint = BlobPathBuilderSpec.buildEndpoint("storageAccount")
    val container = BlobContainerName("storageContainer")
    val refreshedToken = generateTokenExpiration(69L)
    val buffer = Duration.ofMinutes(10L)
    val sasTokenNew = BlobPathBuilderFactorySpec.buildExampleSasToken(refreshedToken)
    val blobTokenGenerator = mock[WSMBlobSasTokenGenerator]
    when(blobTokenGenerator.findBlobSasToken(endpoint, container, buffer)).thenCallRealMethod()
    when(blobTokenGenerator.generateBlobSasToken(endpoint, container)).thenReturn(Success(sasTokenNew))
    when(blobTokenGenerator.getAvailableCachedSasToken(endpoint, container)).thenReturn(None)
    val sas: Try[AzureSasCredential] = blobTokenGenerator.findBlobSasToken(endpoint, container, buffer)
    verify(blobTokenGenerator, times(1)).generateBlobSasToken(endpoint, container)
    verify(blobTokenGenerator, times(1)).putAvailableCachedSasToken(endpoint, container, sasTokenNew)
    verify(blobTokenGenerator, times(1)).getAvailableCachedSasToken(endpoint, container)
    sas shouldBe Success(sasTokenNew)
  }

  it should "Get a sas token from WSM" in {
    val endpoint = BlobPathBuilderSpec.buildEndpoint("lz43a8a3d21540dfd25f5ace")
    val container = BlobContainerName("sc-1337c9a8-854f-4087-b2ba-99558b7171f6")
    val wsmClientProvider = new HttpWorkspaceManagerClientProvider(WorkspaceManagerURL("https://workspace.dsde-dev.broadinstitute.org"))
    val generator = BlobSasTokenGenerator.createBlobTokenGenerator(wsmClientProvider, Some("" +
      "eyJhbGciOiJSUzI1NiIsImtpZCI6Inp4UnJ6aEw4OTBrYi1vVkZOWkdWUFlxSzYtSkI3RFRlTGVvMllTbzRtc1UiLCJ0eXAiOiJKV1QifQ.eyJlbWFpbCI6ImNmcmVpdGFzLjQ1ODNAZ21haWwuY29tIiwiZ2l2ZW5fbmFtZSI6IkNocmlzdGlhbiIsImZhbWlseV9uYW1lIjoiRnJlaXRhcyIsIm5hbWUiOiJDaHJpc3RpYW4gRnJlaXRhcyIsImlkcCI6Imdvb2dsZS5jb20iLCJpZHBfYWNjZXNzX3Rva2VuIjoieWEyOS5hMEFXWTdDa2xvaHJoSS1XaVhaWDFCRGxLc2tqdmJyUGMwWG05Qk4xelEzczBvWDFNTU4yY3NDclhSYXZDaFN4Mkp3OW9lNmN3aFNYdnpsWHpaaDI2aHdBTzV0M3ZMSS1pQjNrX1JFWnJBNGdoSUY5U2hQYWthWkt1SGFOeUtZXzZzN2NvUFhUUmRueFh3YmdHQkFtX3RrQmtoa1E2SUJ0MGFDZ1lLQWIwU0FSTVNGUUcxdERycDdETnl0MktlMVpOOEU4ZWd0UmpUbUEwMTY2IiwiZ29vZ2xlX2lkIjoiMTA0MDU1MTAwMTMzNjg1NDY5Mzg2IiwicGljdHVyZSI6Imh0dHBzOi8vbGgzLmdvb2dsZXVzZXJjb250ZW50LmNvbS9hL0FBY0hUdGNLOGYzQlVqU2Nsb3F2dmN4VU1NU3RzZGdlMWtkSzMxbjQ5ZFc5PXM5Ni1jIiwiZW1haWxfdmVyaWZpZWQiOnRydWUsInN1YiI6ImRhYWY4YzhkLWM4YzQtNGJiMi1hNWMxLTAzOTI3MmU0MWMyOCIsInRpZCI6ImZkMGJjMGVmLTE3NDctNGVlNi1hYjNlLWQ0ZDZiYjg4MmQ0MCIsImF6cCI6ImJiZDA3ZDQzLTAxY2ItNGI2OS04ZmQwLTU3NDZkOWE1YzlmZSIsInZlciI6IjEuMCIsImlhdCI6MTY4NjI2MTg2MywiYXVkIjoiYmJkMDdkNDMtMDFjYi00YjY5LThmZDAtNTc0NmQ5YTVjOWZlIiwiZXhwIjoxNjg2MjY1NDYzLCJpc3MiOiJodHRwczovL3RlcnJhZGV2YjJjLmIyY2xvZ2luLmNvbS9mZDBiYzBlZi0xNzQ3LTRlZTYtYWIzZS1kNGQ2YmI4ODJkNDAvdjIuMC8iLCJuYmYiOjE2ODYyNjE4NjN9.HboOTTWGDT6F8lGbCus0EcSEAQIIJXyYUogBBdispi6-td7TYqxS0xZ3W1sbXrUFFI-rDne0Ck6Cc2GiL4iS-J27If5u47JjmqW3eW2vyX6uylGrnjHWjZWr6t1UDvpD1SHSXk3U24X-zkeEHC5PY_FJ5CHt4t0_a7wFO3fDL9WMmYpsHwiqOlfa2t874O4AspTpA0RvZqw7gdrcxzmRxbUNSJy_cDeOzB4UqgKENfmXZI1azGsSx-0wFfmaPdy_Hyoydf9cXN7KmA9kfBPff9dnj0bdTaH8-LMxj0fHTNRuDnA_LITg9Pxl7f6mA0iSukW_45_XynX6whDxKt2njw" +
      ""), Duration.ofMinutes(10))
    val sas = generator.findBlobSasToken(endpoint, container, Duration.ofMinutes(10))
    sas.isSuccess shouldBe true
  }
}
