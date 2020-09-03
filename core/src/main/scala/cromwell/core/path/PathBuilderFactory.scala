package cromwell.core.path

import akka.actor.ActorSystem
import cromwell.core.{Dispatcher, WorkflowOptions}
import cats.syntax.traverse._

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
      case (a, b) => factories.indexOf(a) < factories.indexOf(b)
    })
    sortedFactories.traverse(_.withOptions(workflowOptions))
  }
}

/**
  * Provide a method that can instantiate a path builder with the specified workflow options.
  */
trait PathBuilderFactory {
  def withOptions(options: WorkflowOptions)(implicit as: ActorSystem, ec: ExecutionContext): Future[PathBuilder]
}
