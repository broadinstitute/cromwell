package cromwell.services.metadata.impl.archiver

import akka.actor.ActorSystem
import com.typesafe.config.Config
import common.Checked
import common.validation.Validation._
import cromwell.core.filesystem.CromwellFileSystems
import cromwell.core.path.PathFactory.PathBuilders
import cromwell.core.path.{PathBuilderFactory, PathFactory}
import cromwell.core.{WorkflowId, WorkflowOptions}
import net.ceedubs.ficus.Ficus._

import scala.concurrent.Await
import scala.concurrent.duration._
import scala.language.postfixOps
import scala.util.Try

final case class ArchiveMetadataConfig(pathBuilders: PathBuilders,
                                       bucket: String,
                                       bucketReadLimit: Int,
                                       interval: FiniteDuration,
                                       debugLogging: Boolean) {
  def makePath(workflowId: WorkflowId)= PathFactory.buildPath(ArchiveMetadataConfig.pathForWorkflow(workflowId, bucket), pathBuilders)
}

object ArchiveMetadataConfig {

  // TODO: Confirm if this makes sense to the users? Should we store /bucket/parent-wf-id/sub-wf-id or /bucket/wf-id ?
  // When deciding keep in mind that workflows can nest to arbitrary depth, mustn't exceed path length limits: /bucket/parent-wf-id/parent-wf-id/parent-wf-id/parent-wf-id/sub-wf-id
  def pathForWorkflow(id: WorkflowId, bucket: String) = s"gs://$bucket/$id/$id.csv"

  def parseConfig(archiveMetadataConfig: Config)(implicit system: ActorSystem): Checked[ArchiveMetadataConfig] = {
    val defaultMaxInterval: FiniteDuration = 5 minutes
    val defaultDebugLogging = true

    for {
      _ <- Try(archiveMetadataConfig.getConfig("filesystems.gcs")).toCheckedWithContext("parse Carboniter 'filesystems.gcs' field from config")
      pathBuilderFactories <- CromwellFileSystems.instance.factoriesFromConfig(archiveMetadataConfig)
      pathBuilders <- Try(Await.result(PathBuilderFactory.instantiatePathBuilders(pathBuilderFactories.values.toList, WorkflowOptions.empty), 60.seconds))
        .toCheckedWithContext("construct Carboniter path builders from factories")
      bucket <- Try(archiveMetadataConfig.getString("bucket")).toCheckedWithContext("parse Carboniter 'bucket' field from config")
      bucketReadLimit <- Try(archiveMetadataConfig.getOrElse[Int]("bucket-read-limit-bytes", 150000000)).toCheckedWithContext("parse Carboniter 'bucket-read-limit-bytes' field from config")
      interval <- Try(archiveMetadataConfig.getOrElse[FiniteDuration]("interval", defaultMaxInterval)).toChecked
      debugLogging <- Try(archiveMetadataConfig.getOrElse("debug-logging", defaultDebugLogging)).toChecked
    } yield ArchiveMetadataConfig(pathBuilders, bucket, bucketReadLimit, interval, debugLogging)
  }
}
