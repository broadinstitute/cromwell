package cromwell.cloudsupport.azure

import com.typesafe.config.{Config}
import org.slf4j.LoggerFactory
import common.exception.MessageAggregation
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

  private val log = LoggerFactory.getLogger("AzureConfiguration")

  final case class AzureConfigurationException(errorMessages: List[String]) extends MessageAggregation {
    override val exceptionContext = "Azure configuration"
  }

  def apply(config: Config): Try[AzureSasCredential] = {
    val azureConfig = config.getConfig("azure")
    val azureSubscription = azureConfig.getString("subscription")
    val blobContainer = azureConfig.getString("container")
    val azureEndpoint = azureConfig.getString("endpoint")

    def parseURI(string: String): Try[URI] = Try(URI.create(UrlEscapers.urlFragmentEscaper().escape(string)))

    def parseStorageAccount(uri: URI): Try[String] = uri.getHost.split("\\.").find(_.nonEmpty)
      .map(Success(_)).getOrElse(Failure(new AzureConfigurationException(List("Could not parse storage account"))))

    val azureProfile = new AzureProfile(AzureEnvironment.AZURE)

    def azureCredentialBuilder = new DefaultAzureCredentialBuilder()
      .authorityHost(azureProfile.getEnvironment.getActiveDirectoryEndpoint)
      .build

    def authenticateWithSubscription(sub: String) = AzureResourceManager.authenticate(azureCredentialBuilder, azureProfile).withSubscription(sub)

    def azure = {log.debug("Authenticating with Azure Subscription"); authenticateWithSubscription(azureSubscription)}

    def findAzureStorageAccount(storageAccountName: String) = azure.storageAccounts.list.asScala.find(_.name.equals(storageAccountName))
      .map(Success(_)).getOrElse(Failure(new Exception("Azure Storage Account not found")))

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

    /**
      * Generate a BlobSasToken by using the local environment azure identity
      * This will use a default subscription if one is not provided.
      *
      * @return an AzureSasCredential for accessing a blob container
      */
    def generateBlobSasToken: Try[AzureSasCredential] = for {
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
    } yield asc

    return generateBlobSasToken
  }

}
