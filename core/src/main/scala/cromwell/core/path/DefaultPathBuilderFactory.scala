package cromwell.core.path

import akka.actor.{ActorRef, ActorSystem}
import cromwell.core.WorkflowOptions
import cromwell.core.path.PathBuilderFactory.PriorityDefault

import scala.concurrent.{ExecutionContext, Future}

case object DefaultPathBuilderFactory extends PathBuilderFactory {
  override def withOptions(options: WorkflowOptions, serviceRegistryActor: ActorRef)(implicit actorSystem: ActorSystem, ec: ExecutionContext) = Future.successful(DefaultPathBuilder)
  val name = "local"
  val tuple = name -> this

  override def priority: Int = PriorityDefault
}
