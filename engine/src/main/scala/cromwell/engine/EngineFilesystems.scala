package cromwell.engine

import akka.actor.ActorSystem
import com.typesafe.config.{Config, ConfigFactory}
import common.validation.Validation._
import cromwell.core.WorkflowOptions
import cromwell.core.filesystem.CromwellFileSystems
import cromwell.core.path.{DefaultPathBuilderFactory, PathBuilder, PathBuilderFactory}
import net.ceedubs.ficus.Ficus._

import scala.concurrent.Future

object EngineFilesystems {
  private val config: Config = ConfigFactory.load

  private val defaultFileSystemFactory: Map[String, PathBuilderFactory] =
    Option(DefaultPathBuilderFactory.tuple)
      .filter(_ => config.as[Boolean]("engine.filesystems.local.enabled"))
      .toMap

  private val pathBuilderFactories: Map[String, PathBuilderFactory] = {
    CromwellFileSystems.instance.factoriesFromConfig(config.as[Config]("engine"))
      .unsafe("Failed to instantiate engine filesystem") ++ defaultFileSystemFactory
  }

  def configuredPathBuilderFactories: List[PathBuilderFactory] = pathBuilderFactories.values.toList

  def pathBuildersForWorkflow(workflowOptions: WorkflowOptions, factories: List[PathBuilderFactory])(implicit as: ActorSystem): Future[List[PathBuilder]] = {
    PathBuilderFactory.instantiatePathBuilders(factories, workflowOptions)
  }
}
