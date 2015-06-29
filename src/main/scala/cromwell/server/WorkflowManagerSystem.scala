package cromwell.server

import akka.actor.ActorSystem
import cromwell.engine.db.DataAccess
import cromwell.engine.workflow.WorkflowManagerActor

trait WorkflowManagerSystem {
  val systemName = "cromwell-system"
  implicit val actorSystem = ActorSystem(systemName)
  lazy val workflowManagerActor = actorSystem.actorOf(WorkflowManagerActor.props(dataAccess))

  def dataAccess: DataAccess
}
