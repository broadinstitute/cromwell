package cromwell.filesystems.sra

import akka.actor.ActorSystem
import com.typesafe.config.Config
import cromwell.core.WorkflowOptions
import cromwell.core.path.PathBuilderFactory

import scala.concurrent.{ExecutionContext, Future}

final case class SraPathBuilderFactory(globalConfig: Config, instanceConfig: Config) extends PathBuilderFactory {
  def withOptions(options: WorkflowOptions)(implicit as: ActorSystem, ec: ExecutionContext): Future[SraPathBuilder] =
    Future.successful(SraPathBuilderFactory.pathBuilder)
}

object SraPathBuilderFactory {
  private lazy val pathBuilder = new SraPathBuilder
}
