package cromwell.filesystems.blob

import akka.actor.ActorSystem
import com.typesafe.config.Config
import cromwell.core.WorkflowOptions
import cromwell.core.path.PathBuilderFactory
import cromwell.core.path.PathBuilderFactory.PriorityBlob

import java.util.UUID
import scala.concurrent.{ExecutionContext, Future}

final case class SubscriptionId(value: UUID) {override def toString: String = value.toString}
final case class BlobContainerName(value: String) {override def toString: String = value}
final case class StorageAccountName(value: String) {override def toString: String = value}
final case class EndpointURL(value: String) {override def toString: String = value}
final case class WorkspaceId(value: UUID) {override def toString: String = value.toString}
final case class ContainerResourceId(value: UUID) {override def toString: String = value.toString}
final case class WorkspaceManagerURL(value: String) {override def toString: String = value}

final case class BlobPathBuilderFactory(globalConfig: Config, instanceConfig: Config, singletonConfig: BlobFileSystemConfigWrapper) extends PathBuilderFactory {

  private val config = singletonConfig.config
  private val container = config.blobContainerName
  private val endpoint = config.endpointURL
  private val subscription = config.subscriptionId
  private val expiryBufferMinutes = config.expiryBufferMinutes

  val blobTokenGenerator: BlobTokenGenerator = config.workspaceManagerConfig.map { wsmConfig =>
    val wsmClient: WorkspaceManagerApiClientProvider = new HttpWorkspaceManagerClientProvider(wsmConfig.url)
    // parameterizing client instead of URL to make injecting mock client possible
    BlobTokenGenerator.createBlobTokenGenerator(container, endpoint, wsmConfig.workspaceId, wsmConfig.containerResourceId, wsmClient, wsmConfig.b2cToken)
  }.getOrElse(
    BlobTokenGenerator.createBlobTokenGenerator(container, endpoint, subscription)
  )

  val fsm: BlobFileSystemManager = BlobFileSystemManager(container, endpoint, expiryBufferMinutes, blobTokenGenerator)

  override def withOptions(options: WorkflowOptions)(implicit as: ActorSystem, ec: ExecutionContext): Future[BlobPathBuilder] = {
    Future {
      new BlobPathBuilder(container, endpoint)(fsm)
    }
  }

  override def priority: Int = PriorityBlob
}
