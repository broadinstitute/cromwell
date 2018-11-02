package cromwell.filesystems.oss

import akka.actor.ActorSystem
import cromwell.core.WorkflowOptions
import cromwell.core.path.PathBuilderFactory
import cromwell.filesystems.oss.nio.OssStorageConfiguration

import scala.concurrent.{ExecutionContext, Future}

final case class OssPathBuilderFactory(configuration: OssStorageConfiguration
                                ) extends PathBuilderFactory {
  def withOptions(options: WorkflowOptions)(implicit as: ActorSystem, ec: ExecutionContext) = {
    Future.successful(OssPathBuilder.fromConfiguration(configuration, options))
  }
}
