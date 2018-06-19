package cromwell.filesystems.gcs

import akka.actor.ActorSystem
import com.typesafe.config.Config
import cromwell.core.WorkflowOptions
import cromwell.core.path.PathBuilderFactory

import scala.concurrent.{ExecutionContext, Future}

case class MockGcsPathBuilderFactory(globalConfig: Config, instanceConfig: Config) extends PathBuilderFactory {
  def withOptions(options: WorkflowOptions)(implicit as: ActorSystem, ec: ExecutionContext) = {
    Future.successful(MockGcsPathBuilder.withOptions(options))
  }
}
