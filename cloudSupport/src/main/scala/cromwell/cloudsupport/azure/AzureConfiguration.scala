package cromwell.cloudsupport.azure

import com.typesafe.config.{Config}
import com.azure.core.credential.AzureSasCredential
import com.azure.core.management.AzureEnvironment
import com.azure.core.management.profile.AzureProfile
import com.azure.identity.DefaultAzureCredentialBuilder
import com.azure.resourcemanager.AzureResourceManager
import com.azure.resourcemanager.storage.models.StorageAccountKey
import com.azure.storage.blob.sas.{BlobContainerSasPermission, BlobServiceSasSignatureValues}
import com.azure.storage.blob.{BlobContainerClient, BlobContainerClientBuilder}
import com.azure.storage.common.StorageSharedKeyCredential
import com.google.common.net.UrlEscapers

import java.net.URI
import java.time.OffsetDateTime
import scala.jdk.CollectionConverters.IterableHasAsScala
import scala.util.{Failure, Success, Try}

final case class AzureConfiguration private (subscription: String, endpoint: String, container: String) {

}

object AzureConfiguration {

  def apply(config: Config): BlobContainerClient = {
    val azureSubscription = config.getString("subscription")
    val blobContainer = config.getString("container")
    val azureEndpoint = config.getString("endpoint")

    def parseURI(string: String): Try[URI] = Try(URI.create(UrlEscapers.urlFragmentEscaper().escape(string)))

    def parseStorageAccount(uri: URI): Try[String] = uri.getHost.split("\\.").find(_.nonEmpty)
      .map(Success(_)).getOrElse(Failure(new Exception("Could not parse storage account")))

    val azureProfile = new AzureProfile(AzureEnvironment.AZURE)

    def azureCredentialBuilder = new DefaultAzureCredentialBuilder()
      .authorityHost(azureProfile.getEnvironment.getActiveDirectoryEndpoint)
      .build

    def authenticateWithSubscription(sub: String) = AzureResourceManager.authenticate(azureCredentialBuilder, azureProfile).withSubscription(sub)

    def azure = authenticateWithSubscription(azureSubscription)

    def findAzureStorageAccount(storageAccountName: String) = azure.storageAccounts.list.asScala.find(_.name.equals(storageAccountName))
      .map(Success(_)).getOrElse(Failure(new Exception("Azure Storage Account not found.")))

    def buildBlobContainerClient(credential: StorageSharedKeyCredential, endpointURL: String, blobContainerName: String): BlobContainerClient = {
      new BlobContainerClientBuilder()
        .credential(credential)
        .endpoint(endpointURL)
        .containerName(blobContainerName)
        .buildClient()
    }

    val bcsp = new BlobContainerSasPermission()
      .setReadPermission(true)
      .setCreatePermission(true)
      .setListPermission(true)
      .setWritePermission(true)

    def generateBlobContainerClient: Try[BlobContainerClient] = for {
      uri <- parseURI(azureEndpoint)
      configuredAccount <- parseStorageAccount(uri)
      azureAccount <- findAzureStorageAccount(configuredAccount)
      keys = azureAccount.getKeys.asScala
      key <- keys.headOption.fold[Try[StorageAccountKey]](Failure(new Exception("Storage account has no keys")))(Success(_))
      first = key.value
      sskc = new StorageSharedKeyCredential(configuredAccount, first)
      bcc = buildBlobContainerClient(sskc, azureEndpoint, blobContainer)
      bsssv = new BlobServiceSasSignatureValues(OffsetDateTime.now.plusDays(1), bcsp)
      asc = new AzureSasCredential(bcc.generateSas(bsssv))
    } yield bcc

    if (generateBlobContainerClient.isFailure) {
      throw new Exception("Failed to generate Blob Container Client.")
    }
    generateBlobContainerClient.get
  }
}
