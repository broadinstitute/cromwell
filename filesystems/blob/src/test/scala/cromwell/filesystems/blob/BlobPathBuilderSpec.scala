package cromwell.filesystems.blob

import com.azure.core.credential.AzureSasCredential
import com.azure.core.management.AzureEnvironment
import com.azure.core.management.profile.AzureProfile
import com.azure.identity.DefaultAzureCredentialBuilder
import com.azure.resourcemanager.AzureResourceManager
import com.azure.storage.blob.BlobContainerClientBuilder
import com.azure.storage.blob.sas.{BlobContainerSasPermission, BlobServiceSasSignatureValues}
import com.azure.storage.common.StorageSharedKeyCredential
import org.scalatest.flatspec.AnyFlatSpec
import org.scalatest.matchers.should.Matchers

import java.nio.file.Files
import java.time.OffsetDateTime
import scala.jdk.CollectionConverters._
import scala.util.{Failure, Success, Try}

object BlobPathBuilderSpec {
  def buildEndpoint(storageAccount: String) = s"https://$storageAccount.blob.core.windows.net"
}

class BlobPathBuilderSpec extends AnyFlatSpec with Matchers{

  it should "parse a URI into a path" in {
    val endpoint = BlobPathBuilderSpec.buildEndpoint("storageAccount")
    val container = "container"
    val evalPath = "/path/to/file"
    val testString = endpoint + "/" + container + evalPath
    BlobPathBuilder.validateBlobPath(testString, container, endpoint) match {
      case BlobPathBuilder.ValidBlobPath(path) => path should equal(evalPath)
      case BlobPathBuilder.UnparsableBlobPath(errorMessage) => fail(errorMessage)
    }
  }

  it should "bad storage account fails causes URI to fail parse into a path" in {
    val endpoint = BlobPathBuilderSpec.buildEndpoint("storageAccount")
    val container = "container"
    val evalPath = "/path/to/file"
    val testString = BlobPathBuilderSpec.buildEndpoint("badStorageAccount") + container + evalPath
    BlobPathBuilder.validateBlobPath(testString, container, endpoint) match {
      case BlobPathBuilder.ValidBlobPath(path) => fail(s"Valid path: $path found when verifying mismatched storage account")
      case BlobPathBuilder.UnparsableBlobPath(errorMessage) => errorMessage.getMessage() should equal(BlobPathBuilder.invalidBlobPathMessage(container, endpoint))
    }
  }

  it should "bad container fails causes URI to fail parse into a path" in {
    val endpoint = BlobPathBuilderSpec.buildEndpoint("storageAccount")
    val container = "container"
    val evalPath = "/path/to/file"
    val testString = endpoint + "badContainer" + evalPath
    BlobPathBuilder.validateBlobPath(testString, container, endpoint) match {
      case BlobPathBuilder.ValidBlobPath(path) => fail(s"Valid path: $path found when verifying mismatched container")
      case BlobPathBuilder.UnparsableBlobPath(errorMessage) => errorMessage.getMessage() should equal(BlobPathBuilder.invalidBlobPathMessage(container, endpoint))
    }
  }

  it should "build a blob path from a test string and read a file" in {
    val storageAccountName = "coaexternalstorage"
    val profile = new AzureProfile(AzureEnvironment.AZURE)
    val azureCredential = new DefaultAzureCredentialBuilder()
      .authorityHost(profile.getEnvironment.getActiveDirectoryEndpoint)
      .build
    val azure = AzureResourceManager.authenticate(azureCredential, profile).withDefaultSubscription()

    val storageAccounts = azure.storageAccounts()
    val storageAccount = storageAccounts
      .list()
      .asScala
      .find(_.name == storageAccountName)

    val storageAccountKeys = storageAccount match {
      case Some(value) => value.getKeys().asScala.map(_.value())
      case _ => fail("Storage Account not found")
    }

    val storageAccountKey = storageAccountKeys.headOption match {
      case Some(value) => value
      case _ => fail("Storage Account has no keys")
    }

    val keyCredential = new StorageSharedKeyCredential(
      storageAccountName,
      storageAccountKey
    )
    val endpoint = BlobPathBuilderSpec.buildEndpoint(storageAccountName)
    val store = "inputs"
    val blobServiceClient = new BlobContainerClientBuilder()
      .credential(keyCredential)
      .endpoint(endpoint)
      .containerName(store)
      .buildClient()

    val blobContainerSasPermission = new BlobContainerSasPermission()
      .setReadPermission(true)
      .setCreatePermission(true)
      .setListPermission(true)
    val blobServiceSasValues = new BlobServiceSasSignatureValues(
      OffsetDateTime.now.plusDays(1),
      blobContainerSasPermission
    )
    val evalPath = "/test/inputFile.txt"
    val sas = blobServiceClient.generateSas(blobServiceSasValues)
    val testString = endpoint + "/" + store + evalPath
    val blobPathTry: Try[BlobPath] = new BlobPathBuilder(new AzureSasCredential(sas), store, endpoint) build testString

    blobPathTry match {
      case Success(blobPath) => {
        blobPath.container should equal(store)
        blobPath.endpoint should equal(endpoint)
        blobPath.pathAsString should equal(testString)
        blobPath.pathWithoutScheme should equal(BlobPathBuilder.parseURI(endpoint).getHost + "/" + store + evalPath)
        val is = Files.newInputStream(blobPath.nioPath)
        val fileText = (is.readAllBytes.map(_.toChar)).mkString
        fileText should include ("This is my test file!!!! Did it work?")
      }
      case Failure(errorMessage) => fail(errorMessage)
    }
  }
}
