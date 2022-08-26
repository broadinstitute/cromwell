package cromwell.filesystems.blob

import akka.actor.ActorSystem
import com.azure.core.credential.AzureSasCredential
import com.azure.core.management.AzureEnvironment
import com.azure.core.management.profile.AzureProfile
import com.azure.identity.DefaultAzureCredentialBuilder
import com.azure.resourcemanager.AzureResourceManager
import com.azure.storage.blob.BlobContainerClientBuilder
import com.azure.storage.blob.nio.AzureFileSystem
import com.azure.storage.blob.sas.{BlobContainerSasPermission, BlobServiceSasSignatureValues}
import com.azure.storage.common.StorageSharedKeyCredential
import com.typesafe.config.Config
import cromwell.core.WorkflowOptions
import cromwell.core.path.PathBuilderFactory
import net.ceedubs.ficus.Ficus._

import java.net.URI
import java.nio.file.{FileSystem, FileSystemNotFoundException, FileSystems}
import java.time.temporal.ChronoUnit
import java.time.{Duration, OffsetDateTime}
import scala.concurrent.{ExecutionContext, Future}
import scala.jdk.CollectionConverters._
import scala.util.{Failure, Try}

final case class BlobFileSystemConfig(config: Config)
final case class BlobPathBuilderFactory(globalConfig: Config, instanceConfig: Config, singletonConfig: BlobFileSystemConfig) extends PathBuilderFactory {
  val sasToken: String = instanceConfig.as[String]("sas-token")
  val container: String = instanceConfig.as[String]("store")
  val endpoint: String = instanceConfig.as[String]("endpoint")
  val workspaceId: Option[String] = instanceConfig.as[Option[String]]("workspace-id")
  val workspaceManagerURL: Option[String] = singletonConfig.config.as[Option[String]]("workspace-manager-url")

  val fsm = FileSystemManager(container, endpoint, 10, workspaceId, workspaceManagerURL)

  override def withOptions(options: WorkflowOptions)(implicit as: ActorSystem, ec: ExecutionContext): Future[BlobPathBuilder] = {
    Future {
      new BlobPathBuilder(fsm, container, endpoint)
    }
  }
}

case class FileSystemManager(container: String,
    endpoint: String,
    preemptionMinutes: Long,
    workspaceId: Option[String] = None,
    workspaceManagerURL: Option[String] = None) {

  var expiry: Option[TokenExpiration] = None
  val blobTokenGenerator: BlobTokenGenerator = BlobTokenGenerator.createBlobTokenGenerator(
    container, endpoint, workspaceId, workspaceManagerURL)

  def buildConfigMap(credential: AzureSasCredential, container: String): Map[String, Object] = {
    Map((AzureFileSystem.AZURE_STORAGE_SAS_TOKEN_CREDENTIAL, credential),
      (AzureFileSystem.AZURE_STORAGE_FILE_STORES, container),
      (AzureFileSystem.AZURE_STORAGE_SKIP_INITIAL_CONTAINER_CHECK, java.lang.Boolean.TRUE))
  }

  def uri = new URI("azb://?endpoint=" + endpoint)

  def retrieveFilesystem(): Try[FileSystem] = {
    synchronized {
      expiry.map(_.hasTokenExpired) match {
        case Some(false) => Try(FileSystems.getFileSystem(uri)) recoverWith {
          // If no filesystem already exists, this will create a new connection, with the provided configs
          case _: FileSystemNotFoundException => blobTokenGenerator.generateAccessToken.flatMap(generateFilesystem(uri, container, _))
        }
        // If the token has expired, OR there is no token record, try to close the FS and regenerate
        case _ => {
            closeFileSystem(uri)
            blobTokenGenerator.generateAccessToken.flatMap(generateFilesystem(uri, container, _))
        }
      }
    }
  }

  def generateFilesystem(uri: URI, container: String, token: AzureSasCredential): Try[FileSystem] = {
    expiry = Some(TokenExpiration(token, Duration.of(preemptionMinutes, ChronoUnit.MINUTES)))
    Try(FileSystems.newFileSystem(uri, buildConfigMap(token, container).asJava))
  }

  def closeFileSystem(uri: URI): Try[Unit] = Try(FileSystems.getFileSystem(uri)).map(_.close)
}

sealed trait BlobTokenGenerator {
  def generateAccessToken: Try[AzureSasCredential]
}

object BlobTokenGenerator {
  def createBlobTokenGenerator(container: String, endpoint: String): BlobTokenGenerator = {
    createBlobTokenGenerator(container, endpoint, None, None)
  }
  def createBlobTokenGenerator(container: String, endpoint: String, workspaceId: Option[String], workspaceManagerURL: Option[String]): BlobTokenGenerator = {
     (container: String, endpoint: String, workspaceId, workspaceManagerURL) match {
       case (container, endpoint, None, None) =>
         NativeBlobTokenGenerator(container, endpoint)
       case (container, endpoint, Some(workspaceId), Some(workspaceManagerURL)) =>
         WSMBlobTokenGenerator(container, endpoint, workspaceId, workspaceManagerURL)
       case _ =>
         throw new Exception("Arguments provided do not match any available BlobTokenGenerator implementation.")
     }
  }
}

case class WSMBlobTokenGenerator(container: String, endpoint: String, workspaceId: String, workspaceManagerURL: String) extends BlobTokenGenerator {
  def generateAccessToken: Try[AzureSasCredential] = Failure(new NotImplementedError)
}

case class NativeBlobTokenGenerator(container: String, endpoint: String) extends BlobTokenGenerator {
  def generateAccessToken: Try[AzureSasCredential] = {
    val storageAccountName = BlobPathBuilder.parseStorageAccount(BlobPathBuilder.parseURI(endpoint)) match {
      case Some(storageAccountName) => storageAccountName
      case _ => throw new Exception("Storage account could not be parsed from endpoint")
    }

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
      case Some(value) => value.getKeys.asScala.map(_.value())
      case _ => throw new Exception("Storage Account not found")
    }

    val storageAccountKey = storageAccountKeys.headOption match {
      case Some(value) => value
      case _ => throw new Exception("Storage Account has no keys")
    }

    val keyCredential = new StorageSharedKeyCredential(
      storageAccountName,
      storageAccountKey
    )
    val blobContainerClient = new BlobContainerClientBuilder()
      .credential(keyCredential)
      .endpoint(endpoint)
      .containerName(container)
      .buildClient()

    val blobContainerSasPermission = new BlobContainerSasPermission()
      .setReadPermission(true)
      .setCreatePermission(true)
      .setListPermission(true)
    val blobServiceSasSignatureValues = new BlobServiceSasSignatureValues(
      OffsetDateTime.now.plusDays(1),
      blobContainerSasPermission
    )

    Try(new AzureSasCredential(blobContainerClient.generateSas(blobServiceSasSignatureValues)))
  }
}
