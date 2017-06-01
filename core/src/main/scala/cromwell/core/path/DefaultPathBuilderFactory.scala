package cromwell.core.path

import akka.actor.ActorSystem
import cromwell.core.WorkflowOptions

import scala.concurrent.{ExecutionContext, Future}

case object DefaultPathBuilderFactory extends PathBuilderFactory {
  override def withOptions(options: WorkflowOptions)(implicit as: ActorSystem, ec: ExecutionContext) = Future.successful(DefaultPathBuilder)
}
