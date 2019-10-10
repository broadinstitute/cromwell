package cromwell.services.metadata.hybridcarbonite

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
import scala.util.Try

final case class HybridCarboniteConfig(enabled: Boolean, pathBuilders: PathBuilders, bucket: String) {
  def makePath(workflowId: WorkflowId)= PathFactory.buildPath(HybridCarboniteConfig.pathForWorkflow(workflowId, bucket), pathBuilders)
}

object HybridCarboniteConfig {

  def pathForWorkflow(id: WorkflowId, bucket: String) = s"gs://$bucket/$id/$id.json"
  
  def parseConfig(carboniterConfig: Config)(implicit system: ActorSystem): Checked[HybridCarboniteConfig] = {
    val enable = carboniterConfig.as[Option[Boolean]]("enabled").getOrElse(false)

    for {
      _ <- Try(carboniterConfig.getConfig("filesystems.gcs")).toCheckedWithContext("parse Carboniter 'filesystems.gcs' field from config")
      pathBuilderFactories <- CromwellFileSystems.instance.factoriesFromConfig(carboniterConfig)
      pathBuilders <- Try(Await.result(PathBuilderFactory.instantiatePathBuilders(pathBuilderFactories.values.toList, WorkflowOptions.empty), 10.seconds))
        .toCheckedWithContext("construct Carboniter path builders from factories")
      bucket <- Try(carboniterConfig.getString("bucket")).toCheckedWithContext("parse Carboniter 'bucket' field from config")
    } yield HybridCarboniteConfig(enable, pathBuilders, bucket)
  }
}
