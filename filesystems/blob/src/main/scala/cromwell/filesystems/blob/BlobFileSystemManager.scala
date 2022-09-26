package cromwell.filesystems.blob

import com.azure.core.credential.AzureSasCredential
import com.azure.core.management.AzureEnvironment
import com.azure.core.management.profile.AzureProfile
import com.azure.identity.DefaultAzureCredentialBuilder
import com.azure.resourcemanager.AzureResourceManager
import com.azure.storage.blob.nio.AzureFileSystem
import com.azure.storage.blob.sas.{BlobContainerSasPermission, BlobServiceSasSignatureValues}
import com.azure.storage.blob.{BlobContainerClient, BlobContainerClientBuilder}
import com.azure.storage.common.StorageSharedKeyCredential

import java.net.URI
import java.nio.file.{FileSystem, FileSystemNotFoundException, FileSystems}
import java.time.temporal.ChronoUnit
import java.time.{Duration, Instant, OffsetDateTime}
import scala.jdk.CollectionConverters._
import scala.util.{Failure, Success, Try}
import com.azure.resourcemanager.storage.models.StorageAccountKey

case class FileSystemAPI() {
  def getFileSystem(uri: URI): Try[FileSystem] = Try(FileSystems.getFileSystem(uri))
  def newFileSystem(uri: URI, config: Map[String, Object]): FileSystem = FileSystems.newFileSystem(uri, config.asJava)
  def closeFileSystem(uri: URI): Option[Unit] = getFileSystem(uri).toOption.map(_.close)
}

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
case class BlobFileSystemManager(
    container: BlobContainerName,
    endpoint: EndpointURL,
    expiryBufferMinutes: Long,
    blobTokenGenerator: BlobTokenGenerator,
    fileSystemAPI: FileSystemAPI = FileSystemAPI(),
    private val initialExpiration: Option[Instant] = None) {
  private var expiry: Option[Instant] = initialExpiration
  val buffer: Duration = Duration.of(expiryBufferMinutes, ChronoUnit.MINUTES)

  def getExpiry: Option[Instant] = expiry
  def uri: URI = BlobFileSystemManager.uri(endpoint)
  def isTokenExpired: Boolean = expiry.exists(BlobFileSystemManager.hasTokenExpired(_, buffer))
  def shouldReopenFilesystem: Boolean = isTokenExpired || expiry.isEmpty
  def retrieveFilesystem(): Try[FileSystem] = {
    synchronized {
      shouldReopenFilesystem match {
        case false => fileSystemAPI.getFileSystem(uri).recoverWith {
          // If no filesystem already exists, this will create a new connection, with the provided configs
          case _: FileSystemNotFoundException => blobTokenGenerator.generateAccessToken.flatMap(generateFilesystem(uri, container, _))
        }
        // If the token has expired, OR there is no token record, try to close the FS and regenerate
        case true =>
          fileSystemAPI.closeFileSystem(uri)
          blobTokenGenerator.generateAccessToken.flatMap(generateFilesystem(uri, container, _))
      }
    }
  }

  private def generateFilesystem(uri: URI, container: BlobContainerName, token: AzureSasCredential): Try[FileSystem] = {
    expiry = BlobFileSystemManager.parseTokenExpiry(token)
    if (expiry.isEmpty) return Failure(new Exception("Could not reopen filesystem, no expiration found"))
    Try(fileSystemAPI.newFileSystem(uri, BlobFileSystemManager.buildConfigMap(token, container)))
  }

}

sealed trait BlobTokenGenerator {def generateAccessToken: Try[AzureSasCredential]}
object BlobTokenGenerator {
  def createBlobTokenGenerator(container: BlobContainerName, endpoint: EndpointURL, subscription: Option[SubscriptionId]): BlobTokenGenerator = {
    createBlobTokenGenerator(container, endpoint, None, None, subscription)
  }
  def createBlobTokenGenerator(container: BlobContainerName,
                               endpoint: EndpointURL,
                               workspaceId: Option[WorkspaceId],
                               workspaceManagerURL: Option[WorkspaceManagerURL],
                               subscription: Option[SubscriptionId]
                              ): BlobTokenGenerator = {
     (container: BlobContainerName, endpoint: EndpointURL, workspaceId, workspaceManagerURL) match {
       case (container, endpoint, None, None) =>
         NativeBlobTokenGenerator(container, endpoint, subscription)
       case (container, endpoint, Some(workspaceId), Some(workspaceManagerURL)) =>
         WSMBlobTokenGenerator(container, endpoint, workspaceId, workspaceManagerURL)
       case _ =>
         throw new Exception("Arguments provided do not match any available BlobTokenGenerator implementation.")
     }
  }
  def createBlobTokenGenerator(container: BlobContainerName, endpoint: EndpointURL): BlobTokenGenerator = createBlobTokenGenerator(container, endpoint, None)
  def createBlobTokenGenerator(container: BlobContainerName, endpoint: EndpointURL, workspaceId: Option[WorkspaceId], workspaceManagerURL: Option[WorkspaceManagerURL]): BlobTokenGenerator =
      createBlobTokenGenerator(container, endpoint, workspaceId, workspaceManagerURL, None)

}

case class WSMBlobTokenGenerator(container: BlobContainerName, endpoint: EndpointURL, workspaceId: WorkspaceId, workspaceManagerURL: WorkspaceManagerURL) extends BlobTokenGenerator {
  def generateAccessToken: Try[AzureSasCredential] = Failure(new NotImplementedError)
}

case class NativeBlobTokenGenerator(container: BlobContainerName, endpoint: EndpointURL, subscription: Option[SubscriptionId] = None) extends BlobTokenGenerator {

  private val azureProfile = new AzureProfile(AzureEnvironment.AZURE)
  private def azureCredentialBuilder = new DefaultAzureCredentialBuilder()
      .authorityHost(azureProfile.getEnvironment.getActiveDirectoryEndpoint)
      .build
  private def authenticateWithSubscription(sub: SubscriptionId) = AzureResourceManager.authenticate(azureCredentialBuilder, azureProfile).withSubscription(sub.value)
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


  def generateAccessToken: Try[AzureSasCredential] = for {
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
