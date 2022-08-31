package cromwell.filesystems.blob

import akka.actor.ActorSystem
import com.azure.core.credential.AzureSasCredential
import com.azure.core.management.AzureEnvironment
import com.azure.core.management.profile.AzureProfile
import com.azure.identity.DefaultAzureCredentialBuilder
import com.azure.resourcemanager.AzureResourceManager
import com.azure.resourcemanager.storage.models.{StorageAccount, StorageAccountKey}
import com.azure.storage.blob.nio.AzureFileSystem
import com.azure.storage.blob.sas.{BlobContainerSasPermission, BlobServiceSasSignatureValues}
import com.azure.storage.blob.{BlobContainerClient, BlobContainerClientBuilder}
import com.azure.storage.common.StorageSharedKeyCredential
import com.typesafe.config.Config
import cromwell.core.WorkflowOptions
import cromwell.core.path.PathBuilderFactory
import net.ceedubs.ficus.Ficus._

import java.net.URI
import java.nio.file.{FileSystem, FileSystemNotFoundException, FileSystems}
import java.time.temporal.ChronoUnit
import java.time.{Duration, Instant, OffsetDateTime}
import scala.concurrent.{ExecutionContext, Future}
import scala.jdk.CollectionConverters._
import scala.util.{Failure, Success, Try}

final case class BlobFileSystemConfig(config: Config)
final case class BlobPathBuilderFactory(globalConfig: Config, instanceConfig: Config, singletonConfig: BlobFileSystemConfig) extends PathBuilderFactory {
  val container: BlobContainerName = BlobContainerName(instanceConfig.as[String]("store"))
  val endpoint: EndpointURL = EndpointURL(instanceConfig.as[String]("endpoint"))
  val workspaceId: Option[WorkspaceId] = instanceConfig.as[Option[String]]("workspace-id").map(WorkspaceId(_))
  val expiryBufferMinutes: Long = instanceConfig.as[Option[Long]]("expiry-buffer-minutes").getOrElse(10)
  val workspaceManagerURL: Option[WorkspaceManagerURL] = singletonConfig.config.as[Option[String]]("workspace-manager-url").map(WorkspaceManagerURL(_))

  val blobTokenGenerator: BlobTokenGenerator = BlobTokenGenerator.createBlobTokenGenerator(
    container, endpoint, workspaceId, workspaceManagerURL)
  val fsm: BlobFileSystemManager = BlobFileSystemManager(container, endpoint, expiryBufferMinutes, blobTokenGenerator)

  override def withOptions(options: WorkflowOptions)(implicit as: ActorSystem, ec: ExecutionContext): Future[BlobPathBuilder] = {
    Future {
      new BlobPathBuilder(container, endpoint)(fsm)
    }
  }
}

final case class BlobContainerName(value: String) {override def toString: String = value}
final case class StorageAccountName(value: String) {override def toString: String = value}
final case class EndpointURL(value: String) {override def toString: String = value}
final case class WorkspaceId(value: String) {override def toString: String = value}
final case class WorkspaceManagerURL(value: String) {override def toString: String = value}

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
case class BlobFileSystemManager(container: BlobContainerName,
    endpoint: EndpointURL,
    expiryBufferMinutes: Long,
    blobTokenGenerator: BlobTokenGenerator,
    fileSystemAPI: FileSystemAPI = FileSystemAPI(),
    initialExpiration: Option[Instant] = None) {
  private var expiry: Option[Instant] = initialExpiration
  val buffer: Duration = Duration.of(expiryBufferMinutes, ChronoUnit.MINUTES)

  def getExpiry: Option[Instant] = expiry
  def uri: URI = BlobFileSystemManager.uri(endpoint)
  def hasTokenExpired: Boolean = expiry.exists(BlobFileSystemManager.hasTokenExpired(_, buffer))
  def retrieveFilesystem(): Try[FileSystem] = {
    synchronized {
      (hasTokenExpired, expiry) match {
        case (false, Some(_)) => fileSystemAPI.getFileSystem(uri) recoverWith {
          // If no filesystem already exists, this will create a new connection, with the provided configs
          case _: FileSystemNotFoundException => blobTokenGenerator.generateAccessToken.flatMap(generateFilesystem(uri, container, _))
        }
        // If the token has expired, OR there is no token record, try to close the FS and regenerate
        case _ =>
          closeFileSystem(uri)
          blobTokenGenerator.generateAccessToken.flatMap(generateFilesystem(uri, container, _))
      }
    }
  }

  private def generateFilesystem(uri: URI, container: BlobContainerName, token: AzureSasCredential): Try[FileSystem] = {
    expiry = BlobFileSystemManager.parseTokenExpiry(token)
    Try(fileSystemAPI.newFileSystem(uri, BlobFileSystemManager.buildConfigMap(token, container)))
  }

  private def closeFileSystem(uri: URI): Try[Unit] = fileSystemAPI.getFileSystem(uri).map(_.close)
}

case class FileSystemAPI() {
  def getFileSystem(uri: URI): Try[FileSystem] = Try(FileSystems.getFileSystem(uri))
  def newFileSystem(uri: URI, config: Map[String, Object]): FileSystem = FileSystems.newFileSystem(uri, config.asJava)
}

sealed trait BlobTokenGenerator {def generateAccessToken: Try[AzureSasCredential]}
object BlobTokenGenerator {
  def createBlobTokenGenerator(container: BlobContainerName, endpoint: EndpointURL): BlobTokenGenerator = {
    createBlobTokenGenerator(container, endpoint, None, None)
  }
  def createBlobTokenGenerator(container: BlobContainerName, endpoint: EndpointURL, workspaceId: Option[WorkspaceId], workspaceManagerURL: Option[WorkspaceManagerURL]): BlobTokenGenerator = {
     (container: BlobContainerName, endpoint: EndpointURL, workspaceId, workspaceManagerURL) match {
       case (container, endpoint, None, None) =>
         NativeBlobTokenGenerator(container, endpoint)
       case (container, endpoint, Some(workspaceId), Some(workspaceManagerURL)) =>
         WSMBlobTokenGenerator(container, endpoint, workspaceId, workspaceManagerURL)
       case _ =>
         throw new Exception("Arguments provided do not match any available BlobTokenGenerator implementation.")
     }
  }
}

case class WSMBlobTokenGenerator(container: BlobContainerName, endpoint: EndpointURL, workspaceId: WorkspaceId, workspaceManagerURL: WorkspaceManagerURL) extends BlobTokenGenerator {
  def generateAccessToken: Try[AzureSasCredential] = Failure(new NotImplementedError)
}

case class NativeBlobTokenGenerator(container: BlobContainerName, endpoint: EndpointURL) extends BlobTokenGenerator {

  private val azureProfile = new AzureProfile(AzureEnvironment.AZURE)
  private def azureCredentialBuilder = new DefaultAzureCredentialBuilder()
      .authorityHost(azureProfile.getEnvironment.getActiveDirectoryEndpoint)
      .build
  private def azure = AzureResourceManager.authenticate(azureCredentialBuilder, azureProfile).withDefaultSubscription()

  private def findAzureStorageAccount(name: StorageAccountName) = azure.storageAccounts.list.asScala.find(_.name.equals(name.value))
      .fold[Try[StorageAccount]](Failure(new Exception("Azure Storage Account not found")))(Success(_))
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


  def generateAccessToken: Try[AzureSasCredential] = for {
    configuredAccount <- BlobPathBuilder.parseStorageAccount(BlobPathBuilder.parseURI(endpoint.value))
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

