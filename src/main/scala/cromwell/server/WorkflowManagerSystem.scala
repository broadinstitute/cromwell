package cromwell.server

import akka.actor.ActorSystem
import cromwell.engine.backend.Backend
import cromwell.engine.workflow.WorkflowManagerActor

trait WorkflowManagerSystem {
  lazy val backend: Backend = WorkflowManagerActor.BackendInstance

  protected def systemName = "cromwell-system"

  protected def newActorSystem(): ActorSystem = ActorSystem(systemName)

  implicit final val actorSystem = newActorSystem()

  // For now there's only one WorkflowManagerActor so no need to dynamically name it
  lazy val workflowManagerActor = actorSystem.actorOf(WorkflowManagerActor.props(backend), "WorkflowManagerActor")
}
