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
import java.nio.file.spi.FileSystemProvider
import java.time.temporal.ChronoUnit
import java.time.{Duration, OffsetDateTime}
import java.util.UUID
import scala.collection.mutable
import scala.jdk.CollectionConverters._
import scala.util.{Failure, Success, Try}

// We encapsulate this functionality here so that we can easily mock it out, to allow for testing without
// actually connecting to Blob storage.
case class AzureFileSystemAPI(private val provider: FileSystemProvider = new AzureFileSystemProvider()) {
  def getFileSystem(uri: URI): Try[AzureFileSystem] = Try(provider.getFileSystem(uri).asInstanceOf[AzureFileSystem])
  def newFileSystem(uri: URI, config: Map[String, Object]): Try[AzureFileSystem] = Try(
    provider.newFileSystem(uri, config.asJava).asInstanceOf[AzureFileSystem]
  )
  def closeFileSystem(uri: URI): Option[Unit] = getFileSystem(uri).toOption.map(_.close)
}

/**
  * The BlobFileSystemManager is an object that is responsible for managing the open filesystem,
  * and refreshing the SAS token that is used to access the blob container containing that filesystem.
  *
  * See BlobSasTokenGenerator for more information on how a SAS token is generated
  */
object BlobFileSystemManager {

  def buildConfigMap(credential: AzureSasCredential, container: BlobContainerName): Map[String, Object] = {
    // Special handling is done here to provide a special key value pair if the placeholder token is provided
    // This is due to the BlobClient requiring an auth token even for public blob paths.
    val sasTuple =
      if (credential == PLACEHOLDER_TOKEN) (AzureFileSystem.AZURE_STORAGE_PUBLIC_ACCESS_CREDENTIAL, PLACEHOLDER_TOKEN)
      else (AzureFileSystem.AZURE_STORAGE_SAS_TOKEN_CREDENTIAL, credential)

    Map(
      sasTuple,
      (AzureFileSystem.AZURE_STORAGE_FILE_STORES, container.value),
      (AzureFileSystem.AZURE_STORAGE_SKIP_INITIAL_CONTAINER_CHECK, java.lang.Boolean.TRUE)
    )
  }
  def combinedEnpointContainerUri(endpoint: EndpointURL, container: BlobContainerName) = new URI(
    "azb://?endpoint=" + endpoint + "/" + container.value
  )

  val PLACEHOLDER_TOKEN = new AzureSasCredential("this-is-a-public-sas")
}

class BlobFileSystemManager(val expiryBufferMinutes: Long,
                            val blobTokenGenerator: BlobSasTokenGenerator,
                            val fileSystemAPI: AzureFileSystemAPI = AzureFileSystemAPI()
) extends LazyLogging {

  def this(config: BlobFileSystemConfig) =
    this(
      config.expiryBufferMinutes,
      BlobSasTokenGenerator.createBlobTokenGeneratorFromConfig(config)
    )

  def this(rawConfig: Config) = this(BlobFileSystemConfig(rawConfig))

  val buffer: Duration = Duration.of(expiryBufferMinutes, ChronoUnit.MINUTES)

  def retrieveFilesystem(endpoint: EndpointURL, container: BlobContainerName): Try[FileSystem] = {
    val uri = BlobFileSystemManager.combinedEnpointContainerUri(endpoint, container)
    synchronized {
      fileSystemAPI.getFileSystem(uri).filter(!_.isExpired(buffer)).recoverWith {
        // If no filesystem already exists, this will create a new connection, with the provided configs
        case _: FileSystemNotFoundException =>
          logger.info(s"Creating new blob filesystem for URI $uri")
          generateFilesystem(uri, container, endpoint)
        case _: NoSuchElementException =>
          // When the filesystem expires, the above filter results in a
          // NoSuchElementException. If expired, close the filesystem
          // and reopen the filesystem with the fresh token
          logger.info(s"Closing & regenerating token for existing blob filesystem at URI $uri")
          fileSystemAPI.closeFileSystem(uri)
          generateFilesystem(uri, container, endpoint)
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
  private def generateFilesystem(uri: URI, container: BlobContainerName, endpoint: EndpointURL): Try[AzureFileSystem] =
    blobTokenGenerator
      .generateBlobSasToken(endpoint, container)
      .flatMap { (token: AzureSasCredential) =>
        fileSystemAPI.newFileSystem(uri, BlobFileSystemManager.buildConfigMap(token, container))
      }
}

sealed trait BlobSasTokenGenerator {
  def generateBlobSasToken(endpoint: EndpointURL, container: BlobContainerName): Try[AzureSasCredential]
}
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
    config.workspaceManagerConfig
      .map { wsmConfig =>
        val wsmClient: WorkspaceManagerApiClientProvider = new HttpWorkspaceManagerClientProvider(wsmConfig.url)

        // WSM-mediated mediated SAS token generator
        // parameterizing client instead of URL to make injecting mock client possible
        BlobSasTokenGenerator.createBlobTokenGenerator(wsmClient, wsmConfig.overrideWsmAuthToken)
      }
      .getOrElse(
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
  def createBlobTokenGenerator(subscription: Option[SubscriptionId]): BlobSasTokenGenerator =
    NativeBlobSasTokenGenerator(subscription)

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
                               overrideWsmAuthToken: Option[String]
  ): BlobSasTokenGenerator =
    new WSMBlobSasTokenGenerator(workspaceManagerClient, overrideWsmAuthToken)

}

case class WSMTerraCoordinates(wsmEndpoint: String, workspaceId: UUID, containerResourceId: UUID)

class WSMBlobSasTokenGenerator(wsmClientProvider: WorkspaceManagerApiClientProvider,
                               overrideWsmAuthToken: Option[String]
) extends BlobSasTokenGenerator {

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
    val wsmAuthToken: Try[String] = getWsmAuth
    container.workspaceId match {
      // If this is a Terra workspace, request a token from WSM
      case Success(workspaceId) =>
        (for {
          wsmAuth <- wsmAuthToken
          wsmAzureResourceClient = wsmClientProvider.getControlledAzureResourceApi(wsmAuth)
          resourceId <- getContainerResourceId(workspaceId, container, Option(wsmAuth))
          sasToken <- wsmAzureResourceClient.createAzureStorageContainerSasToken(workspaceId, resourceId)
        } yield sasToken).recoverWith {
          // If the storage account was still not found in WSM, this may be a public filesystem
          case exception: ApiException if exception.getCode == 404 => Try(BlobFileSystemManager.PLACEHOLDER_TOKEN)
        }
      // Otherwise assume that the container is public and use a placeholder
      // SAS token to bypass the BlobClient authentication requirement
      case Failure(_) => Try(BlobFileSystemManager.PLACEHOLDER_TOKEN)
    }
  }

  private val cachedContainerResourceIds = new mutable.HashMap[BlobContainerName, UUID]()

  // Optionally provide wsmAuth to avoid acquiring it twice in generateBlobSasToken.
  // In the case that the resourceId is not cached and no auth is provided, this function will acquire a new auth as necessary.
  private def getContainerResourceId(workspaceId: UUID,
                                     container: BlobContainerName,
                                     precomputedWsmAuth: Option[String]
  ): Try[UUID] =
    cachedContainerResourceIds.get(container) match {
      case Some(id) => Try(id) // cache hit
      case _ => // cache miss
        val auth: Try[String] = precomputedWsmAuth.map(auth => Try(auth)).getOrElse(getWsmAuth)
        val resourceId = for {
          wsmAuth <- auth
          wsmResourceApi = wsmClientProvider.getResourceApi(wsmAuth)
          resourceId <- wsmResourceApi.findContainerResourceId(workspaceId, container)
        } yield resourceId
        resourceId.map(id => cachedContainerResourceIds.put(container, id)) // NB: Modifying cache state here.
        cachedContainerResourceIds.get(container) match {
          case Some(uuid) => Try(uuid)
          case _ => Failure(new NoSuchElementException("Could not retrieve container resource ID from WSM"))
        }
    }

  private def getWsmAuth: Try[String] =
    overrideWsmAuthToken match {
      case Some(t) => Success(t)
      case None => AzureCredentials.getAccessToken(None).toTry
    }

  private def parseTerraWorkspaceIdFromPath(blobPath: BlobPath): Try[UUID] =
    if (blobPath.container.value.startsWith("sc-")) Try(UUID.fromString(blobPath.container.value.substring(3)))
    else
      Failure(
        new Exception(
          "Could not parse workspace ID from storage container. Are you sure this is a file in a Terra Workspace?"
        )
      )

  /**
   * Return a REST endpoint that will reply with a sas token for the blob storage container associated with the provided blob path.
   * @param blobPath A blob path of a file living in a blob container that WSM knows about (likely a workspace container).
   * @param tokenDuration How long will the token last after being generated. Default is 8 hours. Sas tokens won't last longer than 24h.
   * NOTE: If a blobPath is provided for a file in a container other than what this token generator was constructed for,
   * this function will make two REST requests. Otherwise, the relevant data is already cached locally.
   */
  def getWSMSasFetchEndpoint(blobPath: BlobPath, tokenDuration: Option[Duration] = None): Try[String] = {
    val wsmEndpoint = wsmClientProvider.getBaseWorkspaceManagerUrl
    val lifetimeQueryParameters: String =
      tokenDuration.map(d => s"?sasExpirationDuration=${d.toSeconds.intValue}").getOrElse("")
    val terraInfo: Try[WSMTerraCoordinates] = for {
      workspaceId <- parseTerraWorkspaceIdFromPath(blobPath)
      containerResourceId <- getContainerResourceId(workspaceId, blobPath.container, None)
      coordinates = WSMTerraCoordinates(wsmEndpoint, workspaceId, containerResourceId)
    } yield coordinates
    terraInfo.map { terraCoordinates =>
      s"${terraCoordinates.wsmEndpoint}/api/workspaces/v1/${terraCoordinates.workspaceId.toString}/resources/controlled/azure/storageContainer/${terraCoordinates.containerResourceId.toString}/getSasToken${lifetimeQueryParameters}"
    }
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
  def generateBlobSasToken(endpoint: EndpointURL, container: BlobContainerName): Try[AzureSasCredential] = {
    val c = AzureUtils.buildContainerClientFromLocalEnvironment(container.toString,
                                                                endpoint.toString,
                                                                subscription.map(_.toString)
    )

    c.map { bcc =>
      val bsssv = new BlobServiceSasSignatureValues(OffsetDateTime.now.plusDays(1), bcsp)
      new AzureSasCredential(bcc.generateSas(bsssv))
    }.orElse(Try(BlobFileSystemManager.PLACEHOLDER_TOKEN))
  }
}
