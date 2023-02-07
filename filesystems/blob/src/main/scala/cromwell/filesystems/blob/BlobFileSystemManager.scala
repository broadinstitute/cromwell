package cromwell.filesystems.blob

import com.azure.core.credential.AzureSasCredential
import com.azure.core.management.AzureEnvironment
import com.azure.core.management.profile.AzureProfile
import com.azure.identity.DefaultAzureCredentialBuilder
import com.azure.resourcemanager.AzureResourceManager
import com.azure.resourcemanager.storage.models.StorageAccountKey
import com.azure.storage.blob.nio.AzureFileSystem
import com.azure.storage.blob.sas.{BlobContainerSasPermission, BlobServiceSasSignatureValues}
import com.azure.storage.blob.{BlobContainerClient, BlobContainerClientBuilder}
import com.azure.storage.common.StorageSharedKeyCredential
import com.typesafe.config.Config
import com.typesafe.scalalogging.LazyLogging
import common.validation.Validation._

import java.net.URI
import java.nio.file.{FileSystem, FileSystemNotFoundException, FileSystems}
import java.time.temporal.ChronoUnit
import java.time.{Duration, Instant, OffsetDateTime}
import scala.jdk.CollectionConverters._
import scala.util.{Failure, Success, Try}

// We encapsulate this functionality here so that we can easily mock it out, to allow for testing without
// actually connecting to Blob storage.
case class FileSystemAPI() {
  def getFileSystem(uri: URI): Try[FileSystem] = Try(FileSystems.getFileSystem(uri))
  def newFileSystem(uri: URI, config: Map[String, Object]): FileSystem = FileSystems.newFileSystem(uri, config.asJava)
  def closeFileSystem(uri: URI): Option[Unit] = getFileSystem(uri).toOption.map(_.close)
}
/**
  * The BlobFileSystemManager is an object that is responsible for managing the open filesystem,
  * and refreshing the SAS token that is used to access the blob container containing that filesystem.
  *
  * See BlobSasTokenGenerator for more information on how a SAS token is generated
  */
object BlobFileSystemManager {
  def parseTokenExpiry(token: AzureSasCredential): Option[Instant] = for {
    expiryString <- token.getSignature.split("&").find(_.startsWith("se")).map(_.replaceFirst("se=","")).map(_.replace("%3A", ":"))
    instant = Instant.parse(expiryString)
  } yield instant

  def buildConfigMap(credential: AzureSasCredential, container: BlobContainerName): Map[String, Object] = {
    Map((AzureFileSystem.AZURE_STORAGE_SAS_TOKEN_CREDENTIAL, credential),
      (AzureFileSystem.AZURE_STORAGE_FILE_STORES, container.value),
      (AzureFileSystem.AZURE_STORAGE_SKIP_INITIAL_CONTAINER_CHECK, java.lang.Boolean.TRUE))
  }
  def hasTokenExpired(tokenExpiry: Instant, buffer: Duration): Boolean = Instant.now.plus(buffer).isAfter(tokenExpiry)
  def uri(endpoint: EndpointURL) = new URI("azb://?endpoint=" + endpoint)
}

class BlobFileSystemManager(val container: BlobContainerName,
                            val endpoint: EndpointURL,
                            val expiryBufferMinutes: Long,
                            val blobTokenGenerator: BlobSasTokenGenerator,
                            val fileSystemAPI: FileSystemAPI = FileSystemAPI(),
                            private val initialExpiration: Option[Instant] = None) extends LazyLogging {

  def this(config: BlobFileSystemConfig) = {
    this(
      config.blobContainerName,
      config.endpointURL,
      config.expiryBufferMinutes,
      BlobSasTokenGenerator.createBlobTokenGeneratorFromConfig(config)
    )
  }

  def this(rawConfig: Config) = this(BlobFileSystemConfig(rawConfig))

  val buffer: Duration = Duration.of(expiryBufferMinutes, ChronoUnit.MINUTES)
  private var expiry: Option[Instant] = initialExpiration

  def getExpiry: Option[Instant] = expiry
  def uri: URI = BlobFileSystemManager.uri(endpoint)
  def isTokenExpired: Boolean = expiry.exists(BlobFileSystemManager.hasTokenExpired(_, buffer))
  def shouldReopenFilesystem: Boolean = isTokenExpired || expiry.isEmpty
  def retrieveFilesystem(): Try[FileSystem] = {
    synchronized {
      shouldReopenFilesystem match {
        case false => fileSystemAPI.getFileSystem(uri).recoverWith {
          // If no filesystem already exists, this will create a new connection, with the provided configs
          case _: FileSystemNotFoundException =>
            logger.info(s"Creating new blob filesystem for URI $uri")
            blobTokenGenerator.generateBlobSasToken.flatMap(generateFilesystem(uri, container, _))
        }
        // If the token has expired, OR there is no token record, try to close the FS and regenerate
        case true =>
          logger.info(s"Closing & regenerating token for existing blob filesystem at URI $uri")
          fileSystemAPI.closeFileSystem(uri)
          blobTokenGenerator.generateBlobSasToken.flatMap(generateFilesystem(uri, container, _))
      }
    }
  }

  private def generateFilesystem(uri: URI, container: BlobContainerName, token: AzureSasCredential): Try[FileSystem] = {
    expiry = BlobFileSystemManager.parseTokenExpiry(token)
    if (expiry.isEmpty) return Failure(new Exception("Could not reopen filesystem, no expiration found"))
    Try(fileSystemAPI.newFileSystem(uri, BlobFileSystemManager.buildConfigMap(token, container)))
  }

}

sealed trait BlobSasTokenGenerator { def generateBlobSasToken: Try[AzureSasCredential] }
object BlobSasTokenGenerator {

  /**
    * Creates the correct BlobSasTokenGenerator based on config inputs. This generator is responsible for producing
    * a valid SAS token for use in accessing a Azure blob storage container.
    * Two types of generators can be produced here:
    * > Workspace Manager (WSM) mediated SAS token generator, used to create SAS tokens that allow access for
    * blob containers mediated by the WSM, and is enabled when a WSM config is provided. This is what is intended
    * for use inside Terra
    *    OR
    * > Native SAS token generator, which obtains a valid SAS token from your local environment to reach blob containers
    * your local azure identity has access to and is the default if a WSM config is not found. This is intended for
    * use outside of Terra
    *
    * Both of these generators require an authentication token to authorize the generation of the SAS token.
    * See BlobSasTokenGenerator for more information on how these generators work.
    *
    * @param config A BlobFileSystemConfig object
    * @return An appropriate BlobSasTokenGenerator
    */
  def createBlobTokenGeneratorFromConfig(config: BlobFileSystemConfig): BlobSasTokenGenerator =
    config.workspaceManagerConfig.map { wsmConfig =>
      val wsmClient: WorkspaceManagerApiClientProvider = new HttpWorkspaceManagerClientProvider(wsmConfig.url)

      // WSM-mediated mediated SAS token generator
      // parameterizing client instead of URL to make injecting mock client possible
      BlobSasTokenGenerator.createBlobTokenGenerator(
        config.blobContainerName,
        config.endpointURL,
        wsmConfig.workspaceId,
        wsmConfig.containerResourceId,
        wsmClient,
        wsmConfig.overrideWsmAuthToken
      )
    }.getOrElse(
      // Native SAS token generator
      BlobSasTokenGenerator.createBlobTokenGenerator(config.blobContainerName, config.endpointURL, config.subscriptionId)
    )

  /**
    * Native SAS token generator, uses the DefaultAzureCredentialBuilder in the local environment
    * to produce a SAS token.
    *
    * @param container The BlobContainerName of the blob container to be accessed by the generated SAS token
    * @param endpoint The EndpointURL containing the storage account of the blob container to be accessed by
    * this SAS token
    * @param subscription Optional subscription parameter to use for local authorization.
    * If one is not provided the default subscription is used
    * @return A NativeBlobTokenGenerator, able to produce a valid SAS token for accessing the provided blob
    * container and endpoint locally
    */
  def createBlobTokenGenerator(container: BlobContainerName,
                               endpoint: EndpointURL,
                               subscription: Option[SubscriptionId]): BlobSasTokenGenerator = {
    NativeBlobSasTokenGenerator(container, endpoint, subscription)
  }

  /**
    * WSM-mediated SAS token generator, uses the DefaultAzureCredentialBuilder in the cloud environment
    * to request a SAS token from the WSM to access the given blob container. If an overrideWsmAuthToken
    * is provided this is used instead.
    *
    * @param container The BlobContainerName of the blob container to be accessed by the generated SAS token
    * @param endpoint The EndpointURL containing the storage account of the blob container to be accessed by
    * this SAS token
    * @param workspaceId The WorkspaceId of the account to authenticate against
    * @param containerResourceId The ContainterResourceId of the blob container as WSM knows it
    * @param workspaceManagerClient The client for making requests against WSM
    * @param overrideWsmAuthToken An optional WsmAuthToken used for authenticating against the WSM for a valid
    * SAS token to access the given container and endpoint. This is a dev only option that is only intended
    * for local testing of the WSM interface
    * @return A WSMBlobTokenGenerator, able to produce a valid SAS token for accessing the provided blob
    * container and endpoint that is managed by WSM
    */
  def createBlobTokenGenerator(container: BlobContainerName,
                               endpoint: EndpointURL,
                               workspaceId: WorkspaceId,
                               containerResourceId: ContainerResourceId,
                               workspaceManagerClient: WorkspaceManagerApiClientProvider,
                               overrideWsmAuthToken: Option[String]): BlobSasTokenGenerator = {
    WSMBlobSasTokenGenerator(container, endpoint, workspaceId, containerResourceId, workspaceManagerClient, overrideWsmAuthToken)
  }

}

case class WSMBlobSasTokenGenerator(container: BlobContainerName,
                                    endpoint: EndpointURL,
                                    workspaceId: WorkspaceId,
                                    containerResourceId: ContainerResourceId,
                                    wsmClientProvider: WorkspaceManagerApiClientProvider,
                                    overrideWsmAuthToken: Option[String]) extends BlobSasTokenGenerator {

  /**
    * Generate a BlobSasToken by using the available authorization information
    * If an overrideWsmAuthToken is provided, use this in the wsmClient request
    * Else try to use the environment azure identity to request the SAS token
    *
    * @return an AzureSasCredential for accessing a blob container
    */
  def generateBlobSasToken: Try[AzureSasCredential] = {
    val wsmAuthToken: Try[String] = overrideWsmAuthToken match {
      case Some(t) => Success(t)
      case None => AzureCredentials.getAccessToken(None).toTry
    }

    for {
      wsmAuth <- wsmAuthToken
      wsmClient = wsmClientProvider.getControlledAzureResourceApi(wsmAuth)
      sasToken <- Try(  // Java library throws
        wsmClient.createAzureStorageContainerSasToken(
          workspaceId.value,
          containerResourceId.value,
          null,
          null,
          null,
          null
        ).getToken)
    } yield new AzureSasCredential(sasToken)
  }
}

case class NativeBlobSasTokenGenerator(container: BlobContainerName, endpoint: EndpointURL, subscription: Option[SubscriptionId] = None) extends BlobSasTokenGenerator {
  private val azureProfile = new AzureProfile(AzureEnvironment.AZURE)
  private def azureCredentialBuilder = new DefaultAzureCredentialBuilder()
      .authorityHost(azureProfile.getEnvironment.getActiveDirectoryEndpoint)
      .build
  private def authenticateWithSubscription(sub: SubscriptionId) = AzureResourceManager.authenticate(azureCredentialBuilder, azureProfile).withSubscription(sub.toString)
  private def authenticateWithDefaultSubscription = AzureResourceManager.authenticate(azureCredentialBuilder, azureProfile).withDefaultSubscription()
  private def azure = subscription.map(authenticateWithSubscription(_)).getOrElse(authenticateWithDefaultSubscription)

  private def findAzureStorageAccount(name: StorageAccountName) = azure.storageAccounts.list.asScala.find(_.name.equals(name.value))
      .map(Success(_)).getOrElse(Failure(new Exception("Azure Storage Account not found")))
  private def buildBlobContainerClient(credential: StorageSharedKeyCredential, endpoint: EndpointURL, container: BlobContainerName): BlobContainerClient = {
    new BlobContainerClientBuilder()
        .credential(credential)
        .endpoint(endpoint.value)
        .containerName(container.value)
        .buildClient()
  }
  private val bcsp = new BlobContainerSasPermission()
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
    uri <- BlobPathBuilder.parseURI(endpoint.value)
    configuredAccount <- BlobPathBuilder.parseStorageAccount(uri)
    azureAccount <- findAzureStorageAccount(configuredAccount)
    keys = azureAccount.getKeys.asScala
    key <- keys.headOption.fold[Try[StorageAccountKey]](Failure(new Exception("Storage account has no keys")))(Success(_))
    first = key.value
    sskc = new StorageSharedKeyCredential(configuredAccount.value, first)
    bcc = buildBlobContainerClient(sskc, endpoint, container)
    bsssv = new BlobServiceSasSignatureValues(OffsetDateTime.now.plusDays(1), bcsp)
    asc = new AzureSasCredential(bcc.generateSas(bsssv))
  } yield asc
}
