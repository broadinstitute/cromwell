package cromwell.filesystems.blob

import akka.actor.ActorSystem
import com.typesafe.config.Config
import cromwell.core.WorkflowOptions
import cromwell.core.path.PathBuilderFactory
import cromwell.core.path.PathBuilderFactory.PriorityBlob

import java.util.UUID
import scala.concurrent.{ExecutionContext, Future}
import scala.util.Try

final case class SubscriptionId(value: UUID) { override def toString: String = value.toString }
final case class BlobContainerName(value: String) {
  override def toString: String = value
  lazy val workspaceId: Try[UUID] =
    Try(UUID.fromString(value.replaceFirst("sc-", "")))
}
final case class StorageAccountName(value: String) { override def toString: String = value }
final case class EndpointURL(value: String) {
  override def toString: String = value
  lazy val storageAccountName: Try[StorageAccountName] = {
    val sa = for {
      host <- value.split("//").findLast(_.nonEmpty)
      storageAccountName <- host.split("\\.").find(_.nonEmpty)
    } yield StorageAccountName(storageAccountName)
    sa.toRight(new Exception(s"Storage account name could not be parsed from $value")).toTry
  }
}

final case class WorkspaceManagerURL(value: String) { override def toString: String = value }

final case class BlobPathBuilderFactory(globalConfig: Config, instanceConfig: Config, fsm: BlobFileSystemManager)
    extends PathBuilderFactory {

  override def withOptions(
    options: WorkflowOptions
  )(implicit as: ActorSystem, ec: ExecutionContext): Future[BlobPathBuilder] =
    Future {
      new BlobPathBuilder()(fsm)
    }

  override def priority: Int = PriorityBlob
}
