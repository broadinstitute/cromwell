package cromwell.filesystems.blob

import cromwell.core.path.PathBuilderFactory
import akka.actor.ActorSystem
import cromwell.core.WorkflowOptions
import cromwell.filesystems.blob.BlobPathBuilder
import scala.concurrent.{ExecutionContext, Future}
import com.typesafe.config.Config
import com.azure.core.credential.AzureSasCredential

final case class BlobPathBuilderFactory(globalConfig: Config, instanceConfig: Config) extends PathBuilderFactory {
  override def withOptions(options: WorkflowOptions)(implicit as: ActorSystem, ec: ExecutionContext): Future[BlobPathBuilder] = {
    val sasToken: String = instanceConfig.getString("sasToken")
    val container: String = instanceConfig.getString("store")
    val storageAccount: String = instanceConfig.getString("storageAccount")
    Future {
      new BlobPathBuilder(new AzureSasCredential(sasToken), container, storageAccount)
    }
  }
}
