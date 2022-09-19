package cromwell.filesystems.blob

import akka.actor.ActorSystem
import com.typesafe.config.Config
import cromwell.core.WorkflowOptions
import cromwell.core.path.PathBuilderFactory
import net.ceedubs.ficus.Ficus._

import scala.concurrent.{ExecutionContext, Future}

final case class BlobFileSystemConfig(config: Config)

final case class BlobContainerName(value: String) {override def toString: String = value}
final case class StorageAccountName(value: String) {override def toString: String = value}
final case class EndpointURL(value: String) {override def toString: String = value}
final case class WorkspaceId(value: String) {override def toString: String = value}
final case class WorkspaceManagerURL(value: String) {override def toString: String = value}
final case class BlobPathBuilderFactory(globalConfig: Config, instanceConfig: Config, singletonConfig: BlobFileSystemConfig) extends PathBuilderFactory {
  val subscription: Option[String] = instanceConfig.as[Option[String]]("subscription")
  val container: BlobContainerName = BlobContainerName(instanceConfig.as[String]("container"))
  val endpoint: EndpointURL = EndpointURL(instanceConfig.as[String]("endpoint"))
  val workspaceId: Option[WorkspaceId] = instanceConfig.as[Option[String]]("workspace-id").map(WorkspaceId(_))
  val expiryBufferMinutes: Long = instanceConfig.as[Option[Long]]("expiry-buffer-minutes").getOrElse(10)
  val workspaceManagerURL: Option[WorkspaceManagerURL] = singletonConfig.config.as[Option[String]]("workspace-manager-url").map(WorkspaceManagerURL(_))

  val blobTokenGenerator: BlobTokenGenerator = BlobTokenGenerator.createBlobTokenGenerator(
    container, endpoint, workspaceId, workspaceManagerURL, subscription)
  val fsm: BlobFileSystemManager = BlobFileSystemManager(container, endpoint, expiryBufferMinutes, blobTokenGenerator)

  override def withOptions(options: WorkflowOptions)(implicit as: ActorSystem, ec: ExecutionContext): Future[BlobPathBuilder] = {
    Future {
      new BlobPathBuilder(container, endpoint)(fsm)
    }
  }
}
