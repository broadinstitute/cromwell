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
                                       archiveDelay: FiniteDuration,
                                       instrumentationInterval: FiniteDuration,
                                       debugLogging: Boolean,
                                       batchSize: Long
) {}

object ArchiveMetadataConfig {
  def parseConfig(archiveMetadataConfig: Config)(implicit system: ActorSystem): Checked[ArchiveMetadataConfig] = {
    val defaultMaxInterval: FiniteDuration = 5 minutes
    val defaultArchiveDelay = 365 days
    val defaultInstrumentationInterval = 1 minute
    val defaultDebugLogging = true
    val defaultBatchSize: Long = 1

    for {
      _ <- Try(archiveMetadataConfig.getConfig("filesystems.gcs"))
        .toCheckedWithContext("parse archiver 'filesystems.gcs' field from config")
      pathBuilderFactories <- CromwellFileSystems.instance.factoriesFromConfig(archiveMetadataConfig)
      pathBuilders <- Try(
        Await.result(
          PathBuilderFactory.instantiatePathBuilders(pathBuilderFactories.values.toList, WorkflowOptions.empty),
          60.seconds
        )
      )
        .toCheckedWithContext("construct archiver path builders from factories")
      bucket <- Try(archiveMetadataConfig.getString("bucket"))
        .toCheckedWithContext("parse Carboniter 'bucket' field from config")
      backoffInterval <- Try(
        archiveMetadataConfig.getOrElse[FiniteDuration]("backoff-interval", defaultMaxInterval)
      ).toChecked
      archiveDelay <- Try(archiveMetadataConfig.getOrElse("archive-delay", defaultArchiveDelay)).toChecked
      instrumentationInterval <- Try(
        archiveMetadataConfig.getOrElse("instrumentation-interval", defaultInstrumentationInterval)
      ).toChecked
      debugLogging <- Try(archiveMetadataConfig.getOrElse("debug-logging", defaultDebugLogging)).toChecked
      batchSize <- Try(archiveMetadataConfig.getOrElse("batch-size", defaultBatchSize)).toChecked
    } yield ArchiveMetadataConfig(pathBuilders,
                                  bucket,
                                  backoffInterval,
                                  archiveDelay,
                                  instrumentationInterval,
                                  debugLogging,
                                  batchSize
    )
  }
}
