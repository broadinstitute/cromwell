package cromwell.filesystems.drs

import akka.actor.ActorSystem
import cloud.nio.impl.drs.DrsCloudNioFileSystemProvider
import com.typesafe.config.Config
import cromwell.core.WorkflowOptions
import cromwell.core.path.{PathBuilder, PathBuilderFactory}

import scala.concurrent.{ExecutionContext, Future}


/**
  * Cromwell Wrapper around DrsFileSystems to load the configuration.
  * This class is used as the global configuration class in the drs filesystem
  */
class DrsFileSystemConfig(val config: Config)


class DrsPathBuilderFactory(globalConfig: Config, instanceConfig: Config, singletonConfig: DrsFileSystemConfig) extends PathBuilderFactory {
  override def withOptions(options: WorkflowOptions)(implicit as: ActorSystem, ec: ExecutionContext): Future[PathBuilder] = {
    Future.successful(DrsPathBuilder(new DrsCloudNioFileSystemProvider(singletonConfig.config)))
  }
}
