package cromwell.filesystems.blob

import bio.terra.workspace.client.ApiException
import com.azure.core.credential.AzureSasCredential
import com.azure.storage.blob.nio.{AzureFileSystem, AzureFileSystemProvider}
import com.azure.storage.blob.sas.{BlobContainerSasPermission, BlobServiceSasSignatureValues}
import com.typesafe.config.Config
import com.typesafe.scalalogging.LazyLogging
import common.validation.Validation._
import cromwell.cloudsupport.azure.{AzureCredentials, AzureUtils}

import java.net.URI
import java.nio.file._
import java.time.temporal.ChronoUnit
import java.time.{Duration, Instant, OffsetDateTime}
import java.util.UUID
import scala.jdk.CollectionConverters._
import scala.util.{Failure, Success, Try}

// We encapsulate this functionality here so that we can easily mock it out, to allow for testing without
// actually connecting to Blob storage.
case class FileSystemAPI(private val provider: FileSystemProvider = new AzureFileSystemProvider()) {
  def getFileSystem(uri: URI): Try[FileSystem] = Try(provider.getFileSystem(uri))
  def newFileSystem(uri: URI, config: Map[String, Object]): FileSystem = provider.newFileSystem(uri, config.asJava)
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
    // Special handling is done here to provide a special key value pair if the placeholder token is provided
    // This is due to the BlobClient requiring an auth token even for public blob paths.
    val sasTuple = if (credential == PLACEHOLDER_TOKEN) (AzureFileSystem.AZURE_STORAGE_PUBLIC_ACCESS_CREDENTIAL, PLACEHOLDER_TOKEN)
    else (AzureFileSystem.AZURE_STORAGE_SAS_TOKEN_CREDENTIAL, credential)

    Map(sasTuple, (AzureFileSystem.AZURE_STORAGE_FILE_STORES, container.value),
        (AzureFileSystem.AZURE_STORAGE_SKIP_INITIAL_CONTAINER_CHECK, java.lang.Boolean.TRUE))
  }
  def combinedEnpointContainerUri(endpoint: EndpointURL, container: BlobContainerName) = new URI("azb://?endpoint=" + endpoint + "/" + container.value)

  val PLACEHOLDER_TOKEN = new AzureSasCredential("this-is-a-public-sas")
}

class BlobFileSystemManager(val expiryBufferMinutes: Long,
                            val blobTokenGenerator: BlobSasTokenGenerator,
                            val fileSystemAPI: FileSystemAPI = FileSystemAPI()) extends LazyLogging {

  def this(config: BlobFileSystemConfig) = {
    this(
      config.expiryBufferMinutes,
      BlobSasTokenGenerator.createBlobTokenGeneratorFromConfig(config)
    )
  }

  def this(rawConfig: Config) = this(BlobFileSystemConfig(rawConfig))

  val buffer: Duration = Duration.of(expiryBufferMinutes, ChronoUnit.MINUTES)

  def retrieveFilesystem(endpoint: EndpointURL, container: BlobContainerName): Try[FileSystem] = {
    val uri = BlobFileSystemManager.combinedEnpointContainerUri(endpoint, container)
    synchronized {
      fileSystemAPI.getFileSystem(uri).filter(!_.isExpired(buffer)).recoverWith {
        // If no filesystem already exists, this will create a new connection, with the provided configs
        case _: FileSystemNotFoundException => {
          logger.info(s"Creating new blob filesystem for URI $uri")
          generateFilesystem(uri, container, endpoint)
        }
        case _ : NoSuchElementException => {
          // When the filesystem expires, the above filter results in a
          // NoSuchElementException. If expired, close the filesystem
          // and reopen the filesystem with the fresh token
          logger.info(s"Closing & regenerating token for existing blob filesystem at URI $uri")
          fileSystemAPI.closeFileSystem(uri)
          generateFilesystem(uri, container, endpoint)
        }
      }
    }
  }

  /**
    * Create a new filesystem pointing to a particular container and storage account,
    * generating a SAS token from WSM as needed
    *
    * @param uri a URI formatted to include the scheme, storage account endpoint and container
    * @param container the container to open as a filesystem
    * @param endpoint the endpoint containing the storage account for the container to open
    * @return a try with either the successfully created filesystem, or a failure containing the exception
    */
  private def generateFilesystem(uri: URI, container: BlobContainerName, endpoint: EndpointURL): Try[AzureFileSystem] = {
    blobTokenGenerator.generateBlobSasToken(endpoint, container)
      .flatMap((token: AzureSasCredential) => {
        Try(fileSystemAPI.newFileSystem(uri, BlobFileSystemManager.buildConfigMap(token, container)))
      })
  }
}

sealed trait BlobSasTokenGenerator { def generateBlobSasToken(endpoint: EndpointURL, container: BlobContainerName): Try[AzureSasCredential] }
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
      BlobSasTokenGenerator.createBlobTokenGenerator(wsmClient, wsmConfig.overrideWsmAuthToken)
    }.getOrElse(
      // Native SAS token generator
      BlobSasTokenGenerator.createBlobTokenGenerator(config.subscriptionId)
    )

  /**
    * Native SAS token generator, uses the DefaultAzureCredentialBuilder in the local environment
    * to produce a SAS token.
    *
    * @param subscription Optional subscription parameter to use for local authorization.
    * If one is not provided the default subscription is used
    * @return A NativeBlobTokenGenerator, able to produce a valid SAS token for accessing the provided blob
    * container and endpoint locally
    */
  def createBlobTokenGenerator(subscription: Option[SubscriptionId]): BlobSasTokenGenerator = {
    NativeBlobSasTokenGenerator(subscription)
  }

  /**
    * WSM-mediated SAS token generator, uses the DefaultAzureCredentialBuilder in the cloud environment
    * to request a SAS token from the WSM to access the given blob container. If an overrideWsmAuthToken
    * is provided this is used instead.
    *
    * @param workspaceManagerClient The client for making requests against WSM
    * @param overrideWsmAuthToken An optional WsmAuthToken used for authenticating against the WSM for a valid
    * SAS token to access the given container and endpoint. This is a dev only option that is only intended
    * for local testing of the WSM interface
    * @return A WSMBlobTokenGenerator, able to produce a valid SAS token for accessing the provided blob
    * container and endpoint that is managed by WSM
    */
  def createBlobTokenGenerator(workspaceManagerClient: WorkspaceManagerApiClientProvider,
                               overrideWsmAuthToken: Option[String]): BlobSasTokenGenerator = {
    WSMBlobSasTokenGenerator(workspaceManagerClient, overrideWsmAuthToken)
  }

}

case class WSMBlobSasTokenGenerator(wsmClientProvider: WorkspaceManagerApiClientProvider,
                                    overrideWsmAuthToken: Option[String]) extends BlobSasTokenGenerator {

  /**
    * Generate a BlobSasToken by using the available authorization information
    * If an overrideWsmAuthToken is provided, use this in the wsmClient request
    * Else try to use the environment azure identity to request the SAS token
    * @param endpoint The EndpointURL of the blob container to be accessed by the generated SAS token
    * @param container The BlobContainerName of the blob container to be accessed by the generated SAS token
    *
    * @return an AzureSasCredential for accessing a blob container
    */
  def generateBlobSasToken(endpoint: EndpointURL, container: BlobContainerName): Try[AzureSasCredential] = {
    val wsmAuthToken: Try[String] = overrideWsmAuthToken match {
      case Some(t) => Success(t)
      case None => AzureCredentials.getAccessToken(None).toTry
    }
    container.workspaceId match {
      // If this is a Terra workspace, request a token from WSM
      case Success(workspaceId) => {
        (for {
          wsmAuth <- wsmAuthToken
          wsmAzureResourceClient = wsmClientProvider.getControlledAzureResourceApi(wsmAuth)
          resourceId <- getContainerResourceId(workspaceId, container, wsmAuth)
          sasToken <- wsmAzureResourceClient.createAzureStorageContainerSasToken(workspaceId, resourceId)
        } yield sasToken).recoverWith {
          // If the storage account was still not found in WSM, this may be a public filesystem
          case exception: ApiException if exception.getCode == 404 => Try(BlobFileSystemManager.PLACEHOLDER_TOKEN)
        }
      }
      // Otherwise assume that the container is public and use a placeholder
      // SAS token to bypass the BlobClient authentication requirement
      case Failure(_) => Try(BlobFileSystemManager.PLACEHOLDER_TOKEN)
    }
  }

 def getContainerResourceId(workspaceId: UUID, container: BlobContainerName, wsmAuth : String): Try[UUID] = {
    val wsmResourceClient = wsmClientProvider.getResourceApi(wsmAuth)
    wsmResourceClient.findContainerResourceId(workspaceId, container)
  }
}

case class NativeBlobSasTokenGenerator(subscription: Option[SubscriptionId] = None) extends BlobSasTokenGenerator {
  private val bcsp = new BlobContainerSasPermission()
    .setReadPermission(true)
    .setCreatePermission(true)
    .setListPermission(true)
    .setWritePermission(true)

  /**
    * Generate a BlobSasToken by using the local environment azure identity
    * This will use a default subscription if one is not provided.
    * @param endpoint The EndpointURL of the blob container to be accessed by the generated SAS token
    * @param container The BlobContainerName of the blob container to be accessed by the generated SAS token
    *
    * @return an AzureSasCredential for accessing a blob container
    */
  def generateBlobSasToken(endpoint: EndpointURL, container: BlobContainerName): Try[AzureSasCredential] = for {
    bcc <- AzureUtils.buildContainerClientFromLocalEnvironment(container.toString, endpoint.toString, subscription.map(_.toString))
    bsssv = new BlobServiceSasSignatureValues(OffsetDateTime.now.plusDays(1), bcsp)
    asc = new AzureSasCredential(bcc.generateSas(bsssv))
  } yield asc
}
