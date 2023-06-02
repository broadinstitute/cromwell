package cromwell.filesystems.blob

import com.azure.core.credential.AzureSasCredential
import com.azure.storage.blob.nio.AzureFileSystem
import com.azure.storage.blob.sas.{BlobContainerSasPermission, BlobServiceSasSignatureValues}
import com.typesafe.config.Config
import com.typesafe.scalalogging.LazyLogging
import common.validation.Validation._
import cromwell.cloudsupport.azure.AzureUtils

import java.net.URI
import java.nio.file.{FileSystem, FileSystemNotFoundException, FileSystems}
import java.time.temporal.ChronoUnit
import java.time.{Duration, Instant, OffsetDateTime}
import java.util.concurrent.ConcurrentHashMap
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
  def isSasValid(sas: AzureSasCredential, buffer: Duration): Boolean = parseTokenExpiry(sas).map(!hasTokenExpired(_, buffer)).getOrElse(false)
  def uri(endpoint: EndpointURL) = new URI("azb://?endpoint=" + endpoint)
}

class BlobFileSystemManager(val expiryBufferMinutes: Long,
                            val blobTokenGenerator: BlobSasTokenGenerator,
                            val fileSystemAPI: FileSystemAPI = FileSystemAPI(),
                            private val initialExpiration: Option[Instant] = None,
                            private val initialOpenContainer: Option[String] = None) extends LazyLogging {

  var lastOpenContainer: Option[String] = initialOpenContainer;

  def this(config: BlobFileSystemConfig) = {
    this(
      config.expiryBufferMinutes,
      BlobSasTokenGenerator.createBlobTokenGeneratorFromConfig(config)
    )
  }

  def this(rawConfig: Config) = this(BlobFileSystemConfig(rawConfig))

  val buffer: Duration = Duration.of(expiryBufferMinutes, ChronoUnit.MINUTES)
  private var expiry: Option[Instant] = initialExpiration

  def getExpiry: Option[Instant] = expiry
  def isTokenExpired: Boolean = expiry.exists(BlobFileSystemManager.hasTokenExpired(_, buffer))
  def shouldReopenFilesystem(container: BlobContainerName): Boolean = isTokenExpired || expiry.isEmpty || !lastOpenContainer.contains(container.value)
  def retrieveFilesystem(endpoint: EndpointURL, container: BlobContainerName): Try[FileSystem] = {
    val uri: URI = BlobFileSystemManager.uri(endpoint)
    synchronized {
      shouldReopenFilesystem(container) match {
        case false => fileSystemAPI.getFileSystem(uri).recoverWith {
          // If no filesystem already exists, this will create a new connection, with the provided configs
          case _: FileSystemNotFoundException =>
            logger.info(s"Creating new blob filesystem for URI $uri and container $container, and last container $lastOpenContainer")
            lastOpenContainer = Some(container.value)
            blobTokenGenerator.findBlobSasToken(endpoint, container, buffer).flatMap(generateFilesystem(uri, container, _))
        }
        // If the token has expired, OR there is no token record, try to close the FS and regenerate
        case true =>
          logger.info(s"Closing & regenerating token for existing blob filesystem at URI $uri and container $container, and last container $lastOpenContainer")
          fileSystemAPI.closeFileSystem(uri)
          lastOpenContainer = Some(container.value)
          blobTokenGenerator.findBlobSasToken(endpoint, container, buffer).flatMap(generateFilesystem(uri, container, _))
      }
    }
  }

  private def generateFilesystem(uri: URI, container: BlobContainerName, token: AzureSasCredential): Try[FileSystem] = {
    expiry = BlobFileSystemManager.parseTokenExpiry(token)
    if (expiry.isEmpty) return Failure(new Exception("Could not reopen filesystem, no expiration found"))
    Try(fileSystemAPI.newFileSystem(uri, BlobFileSystemManager.buildConfigMap(token, container)))
  }

}

sealed trait BlobSasTokenGenerator { def findBlobSasToken(endpoint: EndpointURL, container: BlobContainerName, buffer: Duration): Try[AzureSasCredential] }
object BlobSasTokenGenerator {

  /**
    * Creates the correct BlobSasTokenGenerator based on config inputs. This generator is responsible for producing
    * a valid SAS token for use in accessing a Azure blob storage container.
    * Two types of generators can be produced here:
    * > Workspace Manager (WSM) mediated SAS token generator, used to create SAS tokens that allow access for
    * blob containers mediated by the WSM, and is enabled when a WSM config is provided. This SAS token is also cached
    * for future use until it expires. This generator is what is intended for use inside Terra
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
        wsmConfig.workspaceId,
        wsmConfig.containerResourceId,
        wsmClient,
        wsmConfig.overrideWsmAuthToken,
        Duration.ofMinutes(config.expiryBufferMinutes)
      )
    }.getOrElse(
      // Native SAS token generator
      BlobSasTokenGenerator.createBlobTokenGenerator(config.subscriptionId)
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
  def createBlobTokenGenerator(subscription: Option[SubscriptionId]): BlobSasTokenGenerator = {
    NativeBlobSasTokenGenerator(subscription)
  }

  /**
    * WSM-mediated SAS token generator, uses the DefaultAzureCredentialBuilder in the cloud environment
    * to request a SAS token from the WSM to access the given blob container OR retrieve a cached version
    * of a previously fetched SAS token from WSM. If an overrideWsmAuthToken is provided this is used instead.
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
  def createBlobTokenGenerator(workspaceId: WorkspaceId,
                               containerResourceId: ContainerResourceId,
                               workspaceManagerClient: WorkspaceManagerApiClientProvider,
                               overrideWsmAuthToken: Option[String],
                               buffer: Duration): BlobSasTokenGenerator = {
    new WSMBlobSasTokenGenerator(workspaceId, containerResourceId, workspaceManagerClient, overrideWsmAuthToken, buffer)
  }

}

sealed trait SasCacheAvailable;
case class Available(sas: AzureSasCredential) extends SasCacheAvailable;
case class Unavailable() extends SasCacheAvailable;

class WSMBlobSasTokenGenerator(workspaceId: WorkspaceId,
                                    containerResourceId: ContainerResourceId,
                                    wsmClientProvider: WorkspaceManagerApiClientProvider,
                                    overrideWsmAuthToken: Option[String],
                                    buffer: Duration) extends BlobSasTokenGenerator {

  var cachedSasTokens: ConcurrentHashMap[(EndpointURL, BlobContainerName), SasCacheAvailable] = new ConcurrentHashMap[(EndpointURL, BlobContainerName), SasCacheAvailable];
  def getAvailableCachedSasToken(endpoint: EndpointURL, container: BlobContainerName) = cachedSasTokens.getOrDefault((endpoint, container), Unavailable())
  def putAvailableCachedSasToken(endpoint: EndpointURL, container: BlobContainerName, sas: AzureSasCredential) = {
    cachedSasTokens.put((endpoint, container), Available(sas))
  }
  /**
    * Fetch a BlobSasToken from a cache or fall back to using the available authorization information
    * If an overrideWsmAuthToken is provided, use this in the wsmClient request
    * Else try to use the environment azure identity to request the SAS token
    *
    * @return an AzureSasCredential for accessing a blob container
    */
  def findBlobSasToken(endpoint: EndpointURL, container: BlobContainerName, buffer: Duration): Try[AzureSasCredential] = {
     getAvailableCachedSasToken(endpoint, container) match {
      case Available(sas) if BlobFileSystemManager.isSasValid(sas, buffer) => Success(sas)
      // If unavailable or expired refresh SAS cache entry
      case _ => {
        val azureSasTokenTry: Try[AzureSasCredential] = generateBlobSasToken(endpoint, container)
        azureSasTokenTry.toOption.foreach(_ => this.putAvailableCachedSasToken(endpoint, container, _))
        azureSasTokenTry
      }
    }
  }

  /**
    * Requests new SAS token from WSM client for a specified endpoint and container
    * @return
    */
  def generateBlobSasToken(endpoint: EndpointURL, container: BlobContainerName): Try[AzureSasCredential] = {
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
      azureSasToken = new AzureSasCredential(sasToken)
    } yield azureSasToken
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
    *
    * @return an AzureSasCredential for accessing a blob container
    */
  def findBlobSasToken(endpoint: EndpointURL, container: BlobContainerName, buffer: Duration = Duration.ZERO): Try[AzureSasCredential] = for {
    bcc <- AzureUtils.buildContainerClientFromLocalEnvironment(container.toString, endpoint.toString, subscription.map(_.toString))
    bsssv = new BlobServiceSasSignatureValues(OffsetDateTime.now.plusDays(1), bcsp)
    asc = new AzureSasCredential(bcc.generateSas(bsssv))
  } yield asc
}
