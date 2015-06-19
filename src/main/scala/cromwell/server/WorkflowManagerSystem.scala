package cromwell.server

import akka.actor.ActorSystem
import cromwell.engine.db.DataAccess
import cromwell.engine.workflow.WorkflowManagerActor

trait WorkflowManagerSystem {
  val systemName = "cromwell-system"
  def dataAccess: DataAccess
  implicit val actorSystem = ActorSystem(systemName)

  actorSystem.registerOnTermination {actorSystem.log.info(s"$systemName shutting down")}

  val workflowManagerActor = actorSystem.actorOf(WorkflowManagerActor.props(dataAccess))
}
