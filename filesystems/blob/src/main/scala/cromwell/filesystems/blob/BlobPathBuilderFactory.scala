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
  override def withOptions(options: WorkflowOptions)(implicit as: ActorSystem, ec: ExecutionContext): Future[BlobPathBuilder] = {
    val sasToken: String = instanceConfig.getString("sasToken")
    val container: String = instanceConfig.getString("store")
    val endpoint: String = instanceConfig.getString("endpoint")
    Future {
      new BlobPathBuilder(new AzureSasCredential(sasToken), container, endpoint)
    }
  }
}
