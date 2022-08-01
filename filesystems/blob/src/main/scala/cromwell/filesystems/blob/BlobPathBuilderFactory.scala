package cromwell.filesystems.blob

import akka.actor.ActorSystem
import com.azure.core.credential.AzureSasCredential
import com.typesafe.config.Config
import cromwell.core.WorkflowOptions
import cromwell.core.path.PathBuilderFactory
import cromwell.filesystems.blob.BlobPathBuilder
import net.ceedubs.ficus.Ficus._

import scala.concurrent.ExecutionContext
import scala.concurrent.Future

final case class BlobPathBuilderFactory(globalConfig: Config, instanceConfig: Config) extends PathBuilderFactory {
  val sasToken: String = instanceConfig.as[String]("sas-token")
  val container: String = instanceConfig.as[String]("store")
  val endpoint: String = instanceConfig.as[String]("endpoint")
  val workspaceId: String = instanceConfig.as[String]("workspace-id")
  val workspaceManagerURL = globalConfig.as[String]("workspace-manager-url")

  override def withOptions(options: WorkflowOptions)(implicit as: ActorSystem, ec: ExecutionContext): Future[BlobPathBuilder] = {
    Future {
      new BlobPathBuilder(new AzureSasCredential(sasToken), container, endpoint)
    }
  }
}
