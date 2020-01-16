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

final case class HybridCarboniteConfig(enabled: Boolean,
                                       minimumSummaryEntryId: Option[Long],
                                       debugLogging: Boolean,
                                       pathBuilders: PathBuilders,
                                       bucket: String,
                                       freezeScanConfig: HybridCarboniteFreezeScanConfig) {
  def makePath(workflowId: WorkflowId)= PathFactory.buildPath(HybridCarboniteConfig.pathForWorkflow(workflowId, bucket), pathBuilders)
}

final case class HybridCarboniteFreezeScanConfig(initialInterval: FiniteDuration = 5 seconds, maxInterval: FiniteDuration = 5 minutes, multiplier: Double = 1.1)

object HybridCarboniteConfig {

  def pathForWorkflow(id: WorkflowId, bucket: String) = s"gs://$bucket/$id/$id.json"
  
  def parseConfig(carboniterConfig: Config)(implicit system: ActorSystem): Checked[HybridCarboniteConfig] = {
    val enable = carboniterConfig.getOrElse("enabled", false)
    val minimumSummaryEntryIdCheck = carboniterConfig.as[Option[Long]]("minimum-summary-entry-id") match {
      case Some(l) if l < 0 => "Minimum summary entry ID must be above 0. Omit or set to 0 to allow all entries to be summarized.".invalidNelCheck
      case other => other.validNelCheck
    }
    val carboniteDebugLogging = carboniterConfig.getOrElse("debug-logging", true)

    def freezeScanConfig: Checked[HybridCarboniteFreezeScanConfig] = {
      if (carboniterConfig.hasPath("freeze-scan")) {
        val freeze = carboniterConfig.getConfig("freeze-scan")

        val initialInterval = Try { freeze.getOrElse("initial-interval", 5 seconds) } toErrorOr
        val maxInterval = Try { freeze.getOrElse("max-interval", default = 5 minutes) } toErrorOr
        val multiplier = Try { freeze.getOrElse("multiplier", 1.1) } toErrorOr

        val errorOrFreezeScanConfig = (initialInterval, maxInterval, multiplier) mapN {
          case (i, m, x) => HybridCarboniteFreezeScanConfig(i, m, x)
        } contextualizeErrors "parse Carboniter 'freeze-scan' stanza"

        for {
          freezeScanConfig <- errorOrFreezeScanConfig.toEither
          _ <- if (freezeScanConfig.maxInterval > freezeScanConfig.initialInterval) "".validNelCheck else "'max-interval' must be greater than or equal to 'initial-interval' in Carboniter 'freeze-scan' stanza".invalidNelCheck
          _ <- if (freezeScanConfig.multiplier > 1) "".validNelCheck else "'multiplier' must be greater than 1 in Carboniter 'freeze-scan' stanza".invalidNelCheck
        } yield freezeScanConfig
      } else {
        HybridCarboniteFreezeScanConfig().validNelCheck
      }
    }

    for {
      _ <- Try(carboniterConfig.getConfig("filesystems.gcs")).toCheckedWithContext("parse Carboniter 'filesystems.gcs' field from config")
      pathBuilderFactories <- CromwellFileSystems.instance.factoriesFromConfig(carboniterConfig)
      pathBuilders <- Try(Await.result(PathBuilderFactory.instantiatePathBuilders(pathBuilderFactories.values.toList, WorkflowOptions.empty), 60.seconds))
        .toCheckedWithContext("construct Carboniter path builders from factories")
      bucket <- Try(carboniterConfig.getString("bucket")).toCheckedWithContext("parse Carboniter 'bucket' field from config")
      minimumSummaryEntryId <- minimumSummaryEntryIdCheck
      freezeScan <- freezeScanConfig
    } yield HybridCarboniteConfig(enable, minimumSummaryEntryId, carboniteDebugLogging, pathBuilders, bucket, freezeScan)
  }
}
