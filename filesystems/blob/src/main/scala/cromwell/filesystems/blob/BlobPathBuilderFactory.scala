package cromwell.filesystems.blob

import akka.actor.ActorSystem
import com.typesafe.config.Config
import cromwell.core.WorkflowOptions
import cromwell.core.path.PathBuilderFactory
import cromwell.core.path.PathBuilderFactory.PriorityBlob
import net.ceedubs.ficus.Ficus._

import scala.concurrent.{ExecutionContext, Future}

final case class BlobFileSystemConfig(config: Config)

final case class SubscriptionId(value: String) {override def toString: String = value}
final case class BlobContainerName(value: String) {override def toString: String = value}
final case class StorageAccountName(value: String) {override def toString: String = value}
final case class EndpointURL(value: String) {override def toString: String = value}
final case class WorkspaceId(value: String) {override def toString: String = value}
final case class WorkspaceManagerURL(value: String) {override def toString: String = value}
final case class BlobPathBuilderFactory(globalConfig: Config, instanceConfig: Config, singletonConfig: BlobFileSystemConfig) extends PathBuilderFactory {
  val subscription: Option[SubscriptionId] = instanceConfig.as[Option[String]]("subscription").map(SubscriptionId)
  val container: BlobContainerName = BlobContainerName(instanceConfig.as[String]("container"))
  val endpoint: EndpointURL = EndpointURL(instanceConfig.as[String]("endpoint"))
  val workspaceId: Option[WorkspaceId] = instanceConfig.as[Option[String]]("workspace-id").map(WorkspaceId)
  val expiryBufferMinutes: Long = instanceConfig.as[Option[Long]]("expiry-buffer-minutes").getOrElse(10)
  val workspaceManagerURL: Option[WorkspaceManagerURL] = singletonConfig.config.as[Option[String]]("workspace-manager-url").map(WorkspaceManagerURL)
  val b2cToken: Option[String] = Option("asdf") // instanceConfig.as[Option[String]]("b2cToken")

  val blobTokenGenerator: BlobTokenGenerator = (workspaceManagerURL, b2cToken, workspaceId) match {
    case (Some(url), Some(token), Some(workspaceId)) =>
      val wsmClient: WorkspaceManagerApiClientProvider = new HttpWorkspaceManagerClientProvider(url, token)
      // parameterizing client instead of URL to make injecting mock client possible
      BlobTokenGenerator.createBlobTokenGenerator(container, endpoint, workspaceId, wsmClient)
    case _ =>
      BlobTokenGenerator.createBlobTokenGenerator(container, endpoint, subscription)
  }

  val fsm: BlobFileSystemManager = BlobFileSystemManager(container, endpoint, expiryBufferMinutes, blobTokenGenerator)

  override def withOptions(options: WorkflowOptions)(implicit as: ActorSystem, ec: ExecutionContext): Future[BlobPathBuilder] = {
    Future {
      new BlobPathBuilder(container, endpoint)(fsm)
    }
  }

  override def priority: Int = PriorityBlob
}
