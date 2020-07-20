package cromwell.services.metadata.hybridcarbonite

import akka.actor.ActorSystem
import cats.data.Validated._
import cats.syntax.apply._
import cats.syntax.validated._
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
                                       bucketReadLimit: Int,
                                       freezingConfig: MetadataFreezingConfig,
                                       deletionConfig: MetadataDeletionConfig) {
  def makePath(workflowId: WorkflowId)= PathFactory.buildPath(HybridCarboniteConfig.pathForWorkflow(workflowId, bucket), pathBuilders)
}

sealed trait MetadataFreezingConfig
final case class ActiveMetadataFreezingConfig(initialInterval: FiniteDuration,
                                        maxInterval: FiniteDuration,
                                        multiplier: Double,
                                        minimumSummaryEntryId: Option[Long],
                                        debugLogging: Boolean) extends MetadataFreezingConfig
case object InactiveMetadataFreezingConfig extends MetadataFreezingConfig

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
      if (carboniterConfig.hasPath("metadata-freezing")) {
        val defaultInitialInterval: Duration = Duration.Inf
        val defaultMaxInterval = 5 minutes
        val defaultMultiplier = 1.1
        val defaultMinimumSummaryEntryId: Option[Long] = None
        val defaultDebugLogging = true

        val freeze = carboniterConfig.getConfig("metadata-freezing")

        val initialInterval = Try { freeze.getOrElse("initial-interval", defaultInitialInterval) } toErrorOr
        val maxInterval = Try { freeze.getOrElse("max-interval", defaultMaxInterval) } toErrorOr
        val multiplier = Try { freeze.getOrElse("multiplier", defaultMultiplier) } toErrorOr
        val minimumSummaryEntryId = Try { freeze.getOrElse("minimum-summary-entry-id", defaultMinimumSummaryEntryId) } toErrorOr
        val debugLogging = Try { freeze.getOrElse("debug-logging", defaultDebugLogging)}  toErrorOr

        val errorOrFreezingConfig: ErrorOr[MetadataFreezingConfig] = (initialInterval, maxInterval, multiplier, minimumSummaryEntryId, debugLogging) flatMapN {
          case (i, m, x, s, d) =>
            i match {
              case f: FiniteDuration =>
                val summaryCheck = if (s.exists(_ < 0)) "`minimum-summary-entry-id` must be greater than or equal to 0. Omit or set to 0 to allow all entries to be summarized.".invalidNel else "".validNel
                val maxGteInitialCheck = if (f > m) s"'max-interval' $m should be greater than or equal to finite 'initial-interval' $f.".invalidNel else "".validNel
                val multiplierGt1 = if (x > 1) "".validNel else "`multiplier` must be greater than 1.".invalidNel

                (summaryCheck, maxGteInitialCheck, multiplierGt1) mapN {
                  case (_, _, _) => ActiveMetadataFreezingConfig(f, m, x, s, d)
                }
              case _ => InactiveMetadataFreezingConfig.validNel
            }
        } contextualizeErrors "parse Carboniter 'metadata-freezing' stanza"
        errorOrFreezingConfig
      } toEither else InactiveMetadataFreezingConfig.validNelCheck
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
      bucketReadLimit <- Try(carboniterConfig.getOrElse[Int]("bucket-read-limit-bytes", 150000000)).toCheckedWithContext("parse Carboniter 'bucket-read-limit-bytes' field from config")
      freezingConfig <- metadataFreezingConfig
      metadataDeletion <- metadataDeletionConfig
    } yield HybridCarboniteConfig(pathBuilders, bucket, bucketReadLimit, freezingConfig, metadataDeletion)
  }
}
