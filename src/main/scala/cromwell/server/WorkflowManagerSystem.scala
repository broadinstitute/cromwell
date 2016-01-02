package cromwell.server

import akka.actor.ActorSystem
import com.typesafe.config.ConfigFactory
import cromwell.engine.backend.{CromwellBackend, Backend}
import cromwell.engine.workflow.WorkflowManagerActor
import cromwell.engine.backend.BackendType

trait WorkflowManagerSystem {
  protected def systemName = "cromwell-system"
  protected def newActorSystem(): ActorSystem = ActorSystem(systemName)
  implicit final val actorSystem = newActorSystem()
  def backendType: String = ConfigFactory.load.getConfig("backend").getString("backend")
  lazy val backend: Backend = CromwellBackend.initBackend(backendType, actorSystem)
  // For now there's only one WorkflowManagerActor so no need to dynamically name it
  lazy val workflowManagerActor = actorSystem.actorOf(WorkflowManagerActor.props(backend), "WorkflowManagerActor")
}
