package cromwell.core.path

import akka.actor.ActorSystem
import cromwell.core.WorkflowOptions
import cromwell.core.path.PathBuilderFactory.PriorityDefault

import scala.concurrent.{ExecutionContext, Future}

case object DefaultPathBuilderFactory extends PathBuilderFactory {
  override def withOptions(options: WorkflowOptions)(implicit actorSystem: ActorSystem, ec: ExecutionContext) = Future.successful(DefaultPathBuilder)
  val name = "local"
  val tuple = name -> this

  override def priority: Int = PriorityDefault
}
