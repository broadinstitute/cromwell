package cromwell.core.path

import akka.actor.ActorSystem
import cromwell.core.WorkflowOptions

case object DefaultPathBuilderFactory extends PathBuilderFactory {
  override def withOptions(options: WorkflowOptions)(implicit actorSystem: ActorSystem) = DefaultPathBuilder
}
