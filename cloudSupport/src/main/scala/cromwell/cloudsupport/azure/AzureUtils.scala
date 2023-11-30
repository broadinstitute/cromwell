package cromwell.cloudsupport.azure

import com.azure.core.management.AzureEnvironment
import com.azure.core.management.profile.AzureProfile
import com.azure.identity.DefaultAzureCredentialBuilder
import com.azure.resourcemanager.AzureResourceManager
import com.azure.resourcemanager.storage.models.StorageAccountKey
import com.azure.storage.blob.{BlobContainerClient, BlobContainerClientBuilder}
import com.azure.storage.common.StorageSharedKeyCredential
import com.google.common.net.UrlEscapers

import java.net.URI
import scala.jdk.CollectionConverters.IterableHasAsScala
import scala.util.{Failure, Success, Try}

object AzureUtils {

  /**
    * Generates a BlobContainerClient that can interact with the specified container. Authenticates using the local azure client running on the same machine.
    * @param blobContainer     Name of the blob container. Looks something like "my-blob-container".
    * @param azureEndpoint     Azure endpoint of the container. Looks something like https://somedomain.blob.core.windows.net.
    * @param subscription      Azure subscription. A globally unique identifier. If not provided, a default subscription will be used.
    * @return A blob container client capable of interacting with the specified container.
    */
  def buildContainerClientFromLocalEnvironment(blobContainer: String,
                                               azureEndpoint: String,
                                               subscription: Option[String]
  ): Try[BlobContainerClient] = {
    def parseURI(string: String): Try[URI] = Try(URI.create(UrlEscapers.urlFragmentEscaper().escape(string)))
    def parseStorageAccount(uri: URI): Try[String] = uri.getHost
      .split("\\.")
      .find(_.nonEmpty)
      .map(Success(_))
      .getOrElse(Failure(new Exception("Could not parse storage account")))

    val azureProfile = new AzureProfile(AzureEnvironment.AZURE)

    def azureCredentialBuilder = new DefaultAzureCredentialBuilder()
      .authorityHost(azureProfile.getEnvironment.getActiveDirectoryEndpoint)
      .build

    def authenticateWithSubscription(sub: String) =
      AzureResourceManager.authenticate(azureCredentialBuilder, azureProfile).withSubscription(sub)

    def authenticateWithDefaultSubscription =
      AzureResourceManager.authenticate(azureCredentialBuilder, azureProfile).withDefaultSubscription()

    def azure = subscription.map(authenticateWithSubscription(_)).getOrElse(authenticateWithDefaultSubscription)

    def findAzureStorageAccount(storageAccountName: String) = azure.storageAccounts.list.asScala
      .find(_.name.equals(storageAccountName))
      .map(Success(_))
      .getOrElse(Failure(new Exception("Azure Storage Account not found.")))

    def buildBlobContainerClient(credential: StorageSharedKeyCredential,
                                 endpointURL: String,
                                 blobContainerName: String
    ): BlobContainerClient =
      new BlobContainerClientBuilder()
        .credential(credential)
        .endpoint(endpointURL)
        .containerName(blobContainerName)
        .buildClient()

    def generateBlobContainerClient: Try[BlobContainerClient] = for {
      uri <- parseURI(azureEndpoint)
      configuredAccount <- parseStorageAccount(uri)
      azureAccount <- findAzureStorageAccount(configuredAccount)
      keys = azureAccount.getKeys.asScala
      key <- keys.headOption.fold[Try[StorageAccountKey]](Failure(new Exception("Storage account has no keys")))(
        Success(_)
      )
      first = key.value
      sskc = new StorageSharedKeyCredential(configuredAccount, first)
      bcc = buildBlobContainerClient(sskc, azureEndpoint, blobContainer)
    } yield bcc

    generateBlobContainerClient
  }
}
