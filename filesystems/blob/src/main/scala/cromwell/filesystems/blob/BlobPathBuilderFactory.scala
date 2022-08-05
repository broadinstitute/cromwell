package cromwell.filesystems.blob

import akka.actor.ActorSystem
import com.azure.core.management.AzureEnvironment
import com.azure.core.management.profile.AzureProfile
import com.azure.identity.DefaultAzureCredentialBuilder
import com.azure.resourcemanager.AzureResourceManager
import com.azure.storage.blob.BlobServiceClientBuilder
import com.azure.storage.common.StorageSharedKeyCredential
import com.azure.storage.common.sas.{AccountSasPermission, AccountSasResourceType, AccountSasService, AccountSasSignatureValues}
import com.typesafe.config.Config
import cromwell.core.WorkflowOptions
import cromwell.core.path.PathBuilderFactory
import net.ceedubs.ficus.Ficus._

import java.time.OffsetDateTime
import scala.concurrent.ExecutionContext
import scala.concurrent.Future

import scala.jdk.CollectionConverters._

final case class BlobFileSystemConfig(config: Config)
final case class BlobPathBuilderFactory(globalConfig: Config, instanceConfig: Config, singletonConfig: BlobFileSystemConfig) extends PathBuilderFactory {
  val sasToken: String = instanceConfig.as[String]("sas-token")
  val container: String = instanceConfig.as[String]("store")
  val endpoint: String = instanceConfig.as[String]("endpoint")
  val workspaceId: String = instanceConfig.as[String]("workspace-id")
  val workspaceManagerURL: String = singletonConfig.config.as[String]("workspace-manager-url")

  val blobTokenGenerator: BlobTokenGenerator = BlobTokenGenerator.createBlobTokenGenerator(
    container, endpoint, Option(workspaceId), Option(workspaceManagerURL))

  override def withOptions(options: WorkflowOptions)(implicit as: ActorSystem, ec: ExecutionContext): Future[BlobPathBuilder] = {
    Future {
      new BlobPathBuilder(blobTokenGenerator, container, endpoint)
    }
  }
}

sealed trait BlobTokenGenerator {
  def getAccessToken: String
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
         throw new Exception("Bad! >:(")
     }
  }
}

case class WSMBlobTokenGenerator(container: String, endpoint: String, workspaceId: String, workspaceManagerURL: String) extends BlobTokenGenerator {
  def getAccessToken: String = {
    throw new NotImplementedError
  }
}

case class NativeBlobTokenGenerator(container: String, endpoint: String) extends BlobTokenGenerator {
  def getAccessToken: String = {
    val storageAccountName = BlobPathBuilder.parseStorageAccount(BlobPathBuilder.parseURI(endpoint)) match {
      case Some(storageAccountName) => storageAccountName
      case _ => throw new Exception("Storage account could not be parsed from endpoint")
    }

    val profile = new AzureProfile(AzureEnvironment.AZURE)
    val azureCredential = new DefaultAzureCredentialBuilder()
      .authorityHost(profile.getEnvironment.getActiveDirectoryEndpoint)
      .build
    val azure = AzureResourceManager.authenticate(azureCredential, profile).withDefaultSubscription

    val storageAccounts = azure.storageAccounts()
    val storageAccount = storageAccounts
      .list()
      .asScala
      .find(_.name == storageAccountName)

    val storageAccountKeys = storageAccount match {
      case Some(value) => value.getKeys().asScala.map(_.value())
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
    val blobServiceClient = new BlobServiceClientBuilder()
      .credential(keyCredential)
      .endpoint(endpoint)
      .buildClient()

    val accountSasPermission = new AccountSasPermission()
      .setReadPermission(true)
    val services = new AccountSasService()
      .setBlobAccess(true)
    val resourceTypes = new AccountSasResourceType()
      .setObject(true)
      .setContainer(true)
    val accountSasValues = new AccountSasSignatureValues(
      OffsetDateTime.now.plusDays(1),
      accountSasPermission,
      services,
      resourceTypes,
    )

    blobServiceClient.generateAccountSas(accountSasValues)
  }
}
