package cromwell.server

import akka.actor.ActorSystem
import cromwell.engine.backend.Backend
import cromwell.engine.workflow.WorkflowManagerActor

trait WorkflowManagerSystem {
  lazy val backend: Backend = WorkflowManagerActor.BackendInstance

  val systemName = "cromwell-system"
  implicit val actorSystem = ActorSystem(systemName)

  // For now there's only one WorkflowManagerActor so no need to dynamically name it
  lazy val workflowManagerActor = actorSystem.actorOf(WorkflowManagerActor.props(backend), "WorkflowManagerActor")
}
