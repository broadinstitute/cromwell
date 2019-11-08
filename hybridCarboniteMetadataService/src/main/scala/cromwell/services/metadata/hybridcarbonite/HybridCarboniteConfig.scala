package cromwell.services.metadata.hybridcarbonite

import akka.actor.ActorSystem
import com.typesafe.config.Config
import common.Checked
import common.validation.Checked._
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

final case class HybridCarboniteConfig(enabled: Boolean, debugLogging: Boolean, pathBuilders: PathBuilders, bucket: String, freezeScanConfig: HybridCarboniteFreezeScanConfig) {
  def makePath(workflowId: WorkflowId)= PathFactory.buildPath(HybridCarboniteConfig.pathForWorkflow(workflowId, bucket), pathBuilders)
}

final case class HybridCarboniteFreezeScanConfig(initialInterval: FiniteDuration = 5 seconds, maxInterval: FiniteDuration = 5 minutes, multiplier: Double = 1.1)

object HybridCarboniteConfig {

  def pathForWorkflow(id: WorkflowId, bucket: String) = s"gs://$bucket/$id/$id.json"
  
  def parseConfig(carboniterConfig: Config)(implicit system: ActorSystem): Checked[HybridCarboniteConfig] = {
    val enable = carboniterConfig.getOrElse("enabled", false)
    val carboniteDebugLogging = carboniterConfig.getOrElse("debug-logging", true)

    def freezeScanConfig: Checked[HybridCarboniteFreezeScanConfig] = {
      if (carboniterConfig.hasPath("freeze-scan")) {
        val freeze = carboniterConfig.getConfig("freeze-scan")

        val initialInterval = Try { freeze.getOrElse("initial-interval", 5 seconds) }
        val maxInterval = Try { freeze.getOrElse("max-interval", default = 5 minutes) }
        val multiplier = Try { freeze.getOrElse("multiplier", 1.1) }

        import cats.instances.try_._
        import cats.syntax.apply._

        val tryImx = (initialInterval, maxInterval, multiplier) mapN { case (i, m, x) => (i, m, x) }

        for {
          imx <- tryImx.toCheckedWithContext("parse Carboniter 'freeze-scan' stanza")
          (initial, max, multi) = imx
          _ <- if (max > initial) "".validNelCheck else "max-interval must be greater than or equal to initial-interval in Carboniter 'freeze-scan' stanza".invalidNelCheck
          _ <- if (multi > 1) "".validNelCheck else "'multiplier' must be greater than 1 in Carboniter 'freeze-scan' stanza".invalidNelCheck
        } yield HybridCarboniteFreezeScanConfig(initialInterval = initial, maxInterval = max, multiplier = multi)
      } else {
        HybridCarboniteFreezeScanConfig().validNelCheck
      }
    }

    for {
      _ <- Try(carboniterConfig.getConfig("filesystems.gcs")).toCheckedWithContext("parse Carboniter 'filesystems.gcs' field from config")
      pathBuilderFactories <- CromwellFileSystems.instance.factoriesFromConfig(carboniterConfig)
      pathBuilders <- Try(Await.result(PathBuilderFactory.instantiatePathBuilders(pathBuilderFactories.values.toList, WorkflowOptions.empty), 10.seconds))
        .toCheckedWithContext("construct Carboniter path builders from factories")
      bucket <- Try(carboniterConfig.getString("bucket")).toCheckedWithContext("parse Carboniter 'bucket' field from config")
      freezeScan <- freezeScanConfig
    } yield HybridCarboniteConfig(enable, carboniteDebugLogging, pathBuilders, bucket, freezeScan)
  }
}
