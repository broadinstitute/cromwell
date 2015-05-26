package cromwell.server

import akka.actor.ActorSystem
import cromwell.engine.{ActorWorkflowManager, WorkflowManagerActor}

trait WorkflowManagerSystem {
  val systemName = "cromwell-system"
  implicit val actorSystem = ActorSystem(systemName)

  actorSystem.registerOnTermination {actorSystem.log.info(s"$systemName shutting down")}

  val workflowManagerActor = actorSystem.actorOf(ActorWorkflowManager.props)
}
