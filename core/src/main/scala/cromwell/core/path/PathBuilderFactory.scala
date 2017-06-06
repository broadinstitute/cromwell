package cromwell.core.path

import akka.actor.ActorSystem
import cromwell.core.WorkflowOptions

import scala.concurrent.{ExecutionContext, Future}

/**
  * Provide a method that can instantiate a path builder with the specified workflow options.
  */
trait PathBuilderFactory {
  def withOptions(options: WorkflowOptions)(implicit as: ActorSystem, ec: ExecutionContext): Future[PathBuilder]
}
