package cromwell.services.metadata.impl.archiver

import akka.actor.ActorSystem
import com.typesafe.config.Config
import common.Checked
import common.validation.Validation._
import cromwell.core.filesystem.CromwellFileSystems
import cromwell.core.path.PathFactory.PathBuilders
import cromwell.core.path.PathBuilderFactory
import cromwell.core.WorkflowOptions
import net.ceedubs.ficus.Ficus._

import scala.concurrent.Await
import scala.concurrent.duration._
import scala.language.postfixOps
import scala.util.Try

final case class ArchiveMetadataConfig(pathBuilders: PathBuilders,
                                       bucket: String,
                                       backoffInterval: FiniteDuration,
                                       databaseStreamFetchSize: Int,
                                       archiveDelay: FiniteDuration,
                                       debugLogging: Boolean) {
}

object ArchiveMetadataConfig {
  def parseConfig(archiveMetadataConfig: Config)(implicit system: ActorSystem): Checked[ArchiveMetadataConfig] = {
    val defaultMaxInterval: FiniteDuration = 5 minutes
    val defaultArchiveDelay = 365 days
    val defaultDebugLogging = true

    for {
      _ <- Try(archiveMetadataConfig.getConfig("filesystems.gcs")).toCheckedWithContext("parse archiver 'filesystems.gcs' field from config")
      pathBuilderFactories <- CromwellFileSystems.instance.factoriesFromConfig(archiveMetadataConfig)
      pathBuilders <- Try(Await.result(PathBuilderFactory.instantiatePathBuilders(pathBuilderFactories.values.toList, WorkflowOptions.empty), 60.seconds))
        .toCheckedWithContext("construct archiver path builders from factories")
      bucket <- Try(archiveMetadataConfig.getString("bucket")).toCheckedWithContext("parse Carboniter 'bucket' field from config")
      backoffInterval <- Try(archiveMetadataConfig.getOrElse[FiniteDuration]("backoff-interval", defaultMaxInterval)).toChecked
      databaseStreamFetchSize <- Try(archiveMetadataConfig.getOrElse[Int]("database-stream-fetch-size", 100)).toChecked
      archiveDelay <- Try(archiveMetadataConfig.getOrElse("archive-delay", defaultArchiveDelay)).toChecked
      debugLogging <- Try(archiveMetadataConfig.getOrElse("debug-logging", defaultDebugLogging)).toChecked
    } yield ArchiveMetadataConfig(pathBuilders, bucket, backoffInterval, databaseStreamFetchSize, archiveDelay, debugLogging)
  }
}
