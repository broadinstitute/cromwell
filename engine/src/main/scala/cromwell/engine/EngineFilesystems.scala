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

  private val pathBuilderFactories = {
    CromwellFileSystems.instance.factoriesFromConfig(config.as[Config]("engine"))
      .unsafe("Failed to instantiate engine filesystem") ++ defaultFileSystemFactory
  }

  def pathBuildersForWorkflow(workflowOptions: WorkflowOptions)(implicit as: ActorSystem): Future[List[PathBuilder]] = {
    PathBuilderFactory.instantiatePathBuilders(pathBuilderFactories.values.toList, workflowOptions)
  }
}
