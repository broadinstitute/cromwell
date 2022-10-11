package cromwell.core.path

import akka.actor.ActorSystem
import cromwell.core.{Dispatcher, WorkflowOptions}
import cats.syntax.traverse._
import cats.instances.list._
import cats.instances.future._

import scala.concurrent.{ExecutionContext, Future}

object PathBuilderFactory {
  // Given a list of factories, instantiates the corresponding path builders
  def instantiatePathBuilders(factories: List[PathBuilderFactory], workflowOptions: WorkflowOptions)(implicit as: ActorSystem): Future[List[PathBuilder]] = {
    implicit val ec: ExecutionContext = as.dispatchers.lookup(Dispatcher.IoDispatcher)
    // The DefaultPathBuilderFactory always needs to be last.
    // The reason is path builders are tried in order, and the default one is very generous in terms of paths it "thinks" it supports
    // For instance, it will return a Path for a gcs url even though it doesn't really support it
    val sortedFactories = factories.sortWith({
      case (_, DefaultPathBuilderFactory) => true
      case (DefaultPathBuilderFactory, _) => false
      case (a, b) => a.priority < b.priority
    })
    sortedFactories.traverse(_.withOptions(workflowOptions))
  }
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
  def priority: Int = 1000
}
