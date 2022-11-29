package cromwell.filesystems.blob

import akka.actor.ActorSystem
import com.typesafe.config.Config
import cromwell.core.WorkflowOptions
import cromwell.core.path.PathBuilderFactory
import cromwell.core.path.PathBuilderFactory.PriorityBlob
import net.ceedubs.ficus.Ficus._

import scala.concurrent.{ExecutionContext, Future}

final case class BlobFileSystemConfig(config: Config)

// WSM config is needed for accessing WSM-managed blob containers created in Terra workspaces.
// If the identity executing Cromwell has native access to the blob container, this can be ignored.
final case class WorkspaceManagerConfig(
    url: WorkspaceManagerURL,
    workspaceId: WorkspaceId,
    containerResourceId: ContainerResourceId,
    b2cToken: Option[String] // dev-only
)

object WorkspaceManagerConfig {
  def apply(wsmConfig: Config): WorkspaceManagerConfig =
    WorkspaceManagerConfig(
      WorkspaceManagerURL(wsmConfig.as[String]("url")),
      WorkspaceId(wsmConfig.as[String]("workspace-id")),
      ContainerResourceId(wsmConfig.as[String]("container-resource-id")),
      wsmConfig.as[Option[String]]("b2cToken")
    )
}

final case class SubscriptionId(value: String) {override def toString: String = value}
final case class BlobContainerName(value: String) {override def toString: String = value}
final case class StorageAccountName(value: String) {override def toString: String = value}
final case class EndpointURL(value: String) {override def toString: String = value}
final case class WorkspaceId(value: String) {override def toString: String = value}
final case class ContainerResourceId(value: String) {override def toString: String = value}
final case class WorkspaceManagerURL(value: String) {override def toString: String = value}

final case class BlobPathBuilderFactory(globalConfig: Config, instanceConfig: Config, singletonConfig: BlobFileSystemConfig) extends PathBuilderFactory {
  val subscription: Option[SubscriptionId] = singletonConfig.config.as[Option[String]]("subscription").map(SubscriptionId)
  val container: BlobContainerName = BlobContainerName(singletonConfig.config.as[String]("container"))
  val endpoint: EndpointURL = EndpointURL(singletonConfig.config.as[String]("endpoint"))
  val expiryBufferMinutes: Long = singletonConfig.config.as[Option[Long]]("expiry-buffer-minutes").getOrElse(10)
  val workspaceManagerConfig: Option[WorkspaceManagerConfig] =
    if (singletonConfig.config.hasPath("workspace-manager"))
      Option(WorkspaceManagerConfig(singletonConfig.config.getConfig("workspace-manager")))
    else None

  val blobTokenGenerator: BlobTokenGenerator = workspaceManagerConfig.map { wsmConfig =>
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
