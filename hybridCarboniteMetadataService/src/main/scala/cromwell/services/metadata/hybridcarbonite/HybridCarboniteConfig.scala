package cromwell.services.metadata.hybridcarbonite

import akka.actor.ActorSystem
import com.typesafe.config.Config
import common.validation.Validation._
import cromwell.core.{WorkflowId, WorkflowOptions}
import cromwell.core.filesystem.CromwellFileSystems
import cromwell.core.path.{PathBuilderFactory, PathFactory}
import cromwell.core.path.PathFactory.PathBuilders

import scala.concurrent.Await
import scala.concurrent.duration._

final case class HybridCarboniteConfig(pathBuilders: PathBuilders, bucket: String) {
  def makePath(workflowId: WorkflowId)= PathFactory.buildPath(HybridCarboniteConfig.pathForWorkflow(workflowId, bucket), pathBuilders)
}

object HybridCarboniteConfig {

  def pathForWorkflow(id: WorkflowId, bucket: String) = s"gs://$bucket/$id/$id.json"

  // TODO: Make this an IO/Future/??? based method vs this unsafe monstrosity?
  def make(carboniterConfig: Config)(implicit system: ActorSystem): HybridCarboniteConfig = {
    lazy val pathBuilderFactories = CromwellFileSystems.instance.factoriesFromConfig(carboniterConfig).unsafe("Failed to instantiate engine filesystem")

    lazy val pathBuilders = {
      Await.result(PathBuilderFactory.instantiatePathBuilders(pathBuilderFactories.values.toList, WorkflowOptions.empty), 10.seconds)
    }

    lazy val bucket = carboniterConfig.getString("bucket")

    HybridCarboniteConfig(pathBuilders, bucket)
  }
}
