package cromwell.services.metadata.hybridcarbonite

import akka.actor.ActorSystem
import cats.syntax.apply._
import com.typesafe.config.Config
import common.Checked
import common.validation.Checked._
import common.validation.ErrorOr._
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

final case class HybridCarboniteConfig(pathBuilders: PathBuilders,
                                       bucket: String,
                                       freezingConfig: MetadataFreezingConfig,
                                       deletionConfig: MetadataDeletionConfig) {
  def makePath(workflowId: WorkflowId)= PathFactory.buildPath(HybridCarboniteConfig.pathForWorkflow(workflowId, bucket), pathBuilders)
}

final case class MetadataFreezingConfig(initialInterval: Duration,
                                        maxInterval: FiniteDuration,
                                        multiplier: Double,
                                        minimumSummaryEntryId: Option[Long],
                                        debugLogging: Boolean) {
  lazy val enabled = initialInterval.isFinite()
}

sealed trait MetadataDeletionConfig
final case class ActiveMetadataDeletionConfig(initialDelay: FiniteDuration,
                                              interval: FiniteDuration,
                                              batchSize: Long = 200L,
                                              delayAfterWorkflowCompletion: FiniteDuration = 24 hours) extends MetadataDeletionConfig
case object InactiveMetadataDeletionConfig extends MetadataDeletionConfig

object HybridCarboniteConfig {

  def pathForWorkflow(id: WorkflowId, bucket: String) = s"gs://$bucket/$id/$id.json"

  def parseConfig(carboniterConfig: Config)(implicit system: ActorSystem): Checked[HybridCarboniteConfig] = {

    def metadataFreezingConfig: Checked[MetadataFreezingConfig] = {
      val defaultInitialInterval: Duration = Duration.Inf
      val defaultMaxInterval = 5 minutes
      val defaultMultiplier = 1.1
      val defaultMinimumSummaryEntryId: Option[Long] = None
      val defaultDebugLogging = true
      if (carboniterConfig.hasPath("metadata-freezing")) {
        val freeze = carboniterConfig.getConfig("metadata-freezing")

        val initialInterval = Try { freeze.getOrElse("initial-interval", defaultInitialInterval) } toErrorOr
        val maxInterval = Try { freeze.getOrElse("max-interval", defaultMaxInterval) } toErrorOr
        val multiplier = Try { freeze.getOrElse("multiplier", defaultMultiplier) } toErrorOr
        val minimumSummaryEntryId = Try { freeze.getOrElse("minimum-summary-entry-id", defaultMinimumSummaryEntryId) } toErrorOr
        val debugLogging = Try { freeze.getOrElse("debug-logging", defaultDebugLogging)}  toErrorOr

        val errorOrFreezingConfig = (initialInterval, maxInterval, multiplier, minimumSummaryEntryId, debugLogging) mapN {
          case (i, m, x, s, d) => MetadataFreezingConfig(i, m, x, s, d)
        } contextualizeErrors "parse Carboniter 'metadata-freezing' stanza"

        for {
          freezingConfig <- errorOrFreezingConfig.toEither
          _ <- if (freezingConfig.minimumSummaryEntryId.exists(_ < 0)) "`metadata-freezing.minimum-summary-entry-id` must be greater than or equal to 0. Omit or set to 0 to allow all entries to be summarized.".invalidNelCheck else "".validNelCheck
          _ <- if (freezingConfig.initialInterval.isFinite() && freezingConfig.maxInterval > freezingConfig.initialInterval) "".validNelCheck else "'max-interval' must be greater than or equal to a finite 'initial-interval' in Carboniter 'metadata-freezing' stanza".invalidNelCheck
          _ <- if (freezingConfig.multiplier > 1) "".validNelCheck else "`metadata-freezing.multiplier` must be greater than 1 in Carboniter 'metadata-freezing' stanza".invalidNelCheck
        } yield freezingConfig
      } else {
        MetadataFreezingConfig(
          defaultInitialInterval,
          defaultMaxInterval,
          defaultMultiplier,
          defaultMinimumSummaryEntryId,
          defaultDebugLogging).validNelCheck
      }
    }

    def metadataDeletionConfig: Checked[MetadataDeletionConfig] = {
      if (carboniterConfig.hasPath("metadata-deletion")) {
        val metadataDeletion = carboniterConfig.getConfig("metadata-deletion")

        val initialDelay = Try { metadataDeletion.getOrElse[Duration]("initial-delay", Duration.Inf) } toErrorOr
        val interval = Try { metadataDeletion.getOrElse[Duration]("interval", Duration.Inf) } toErrorOr
        val batchSize = Try { metadataDeletion.getOrElse("batchSize", default = 200L) } toErrorOr
        val delayAfterWorkflowCompletion = Try { metadataDeletion.getOrElse("delay-after-workflow-completion", 24 hours) } toErrorOr

        val errorOrMetadataDeletionConfig: ErrorOr[MetadataDeletionConfig] = (initialDelay, interval, batchSize, delayAfterWorkflowCompletion) mapN { case (id, i, b, d) =>
          if (id.isFinite() && i.isFinite()) {
            ActiveMetadataDeletionConfig(
              id.asInstanceOf[FiniteDuration],
              i.asInstanceOf[FiniteDuration],
              b,
              d
            )
          } else InactiveMetadataDeletionConfig
        } contextualizeErrors "parse Carboniter 'metadata-deletion' stanza"
        errorOrMetadataDeletionConfig.toEither
      } else {
        InactiveMetadataDeletionConfig.validNelCheck
      }
    }

    for {
      _ <- Try(carboniterConfig.getConfig("filesystems.gcs")).toCheckedWithContext("parse Carboniter 'filesystems.gcs' field from config")
      pathBuilderFactories <- CromwellFileSystems.instance.factoriesFromConfig(carboniterConfig)
      pathBuilders <- Try(Await.result(PathBuilderFactory.instantiatePathBuilders(pathBuilderFactories.values.toList, WorkflowOptions.empty), 60.seconds))
        .toCheckedWithContext("construct Carboniter path builders from factories")
      bucket <- Try(carboniterConfig.getString("bucket")).toCheckedWithContext("parse Carboniter 'bucket' field from config")
      freezingConfig <- metadataFreezingConfig
      metadataDeletion <- metadataDeletionConfig
    } yield HybridCarboniteConfig(pathBuilders, bucket, freezingConfig, metadataDeletion)
  }
}
