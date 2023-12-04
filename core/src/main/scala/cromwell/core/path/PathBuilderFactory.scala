package cromwell.core.path

import akka.actor.ActorSystem
import cromwell.core.{Dispatcher, WorkflowOptions}
import cats.syntax.traverse._
import cats.instances.list._
import cats.instances.future._
import cromwell.core.path.PathBuilderFactory.PriorityStandard

import scala.concurrent.{ExecutionContext, Future}

object PathBuilderFactory {
  // Given a list of factories, instantiates the corresponding path builders
  def instantiatePathBuilders(factories: List[PathBuilderFactory], workflowOptions: WorkflowOptions)(implicit
    as: ActorSystem
  ): Future[List[PathBuilder]] = {
    implicit val ec: ExecutionContext = as.dispatchers.lookup(Dispatcher.IoDispatcher)
    val sortedFactories = factories.sortBy(_.priority)
    sortedFactories.traverse(_.withOptions(workflowOptions))
  }

  val PriorityBlob =
    100 // High priority to evaluate first, because blob files may inadvertently match other filesystems
  val PriorityStandard = 1000
  val PriorityDefault = 10000 // "Default" is a fallback, evaluate last
}

/**
  * Provide a method that can instantiate a path builder with the specified workflow options.
  */
trait PathBuilderFactory {
  def withOptions(options: WorkflowOptions)(implicit as: ActorSystem, ec: ExecutionContext): Future[PathBuilder]

  /**
    * Candidate filesystems are considered in a stable order, as some requests may match multiple filesystems.
    * To customize this order, the priority of a filesystem may be adjusted. Lower number == higher priority.
    * @return This filesystem's priority
    */
  def priority: Int = PriorityStandard
}
