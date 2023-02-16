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

final case class BlobPathBuilderFactory(globalConfig: Config, instanceConfig: Config, fsm: BlobFileSystemManager) extends PathBuilderFactory {

  override def withOptions(options: WorkflowOptions)(implicit as: ActorSystem, ec: ExecutionContext): Future[BlobPathBuilder] = {
    Future {
      new BlobPathBuilder(fsm.container, fsm.endpoint)(fsm)
    }
  }

  override def priority: Int = PriorityBlob
}
