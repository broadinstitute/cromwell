package cromwell.filesystems.blob

import akka.actor.ActorSystem
import com.azure.core.credential.AzureSasCredential
import com.typesafe.config.Config
import cromwell.core.WorkflowOptions
import cromwell.core.path.PathBuilderFactory
import cromwell.filesystems.blob.BlobPathBuilder

import scala.concurrent.ExecutionContext
import scala.concurrent.Future

final case class BlobPathBuilderFactory(globalConfig: Config, instanceConfig: Config) extends PathBuilderFactory {
  val sasToken: String = instanceConfig.getString("filesystems.blob.instance.sas-token")
  val container: String = instanceConfig.getString("filesystems.blob.instance.store")
  val endpoint: String = instanceConfig.getString("filesystems.blob.instance.endpoint")
  val workspaceId: String = instanceConfig.getString("filesystems.blob.instance.workspaceId")
  val workspaceManagerURL = globalConfig.getString("filesystems.blob.global.workspace-manager-url")

  override def withOptions(options: WorkflowOptions)(implicit as: ActorSystem, ec: ExecutionContext): Future[BlobPathBuilder] = {
    Future {
      new BlobPathBuilder(new AzureSasCredential(sasToken), container, endpoint)
    }
  }
}
