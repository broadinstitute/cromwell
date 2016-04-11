
package cromwell.server

import akka.actor.ActorSystem
import cromwell.engine.backend.{BackendConfiguration, CromwellBackend}
import cromwell.engine.workflow.WorkflowManagerActor

trait WorkflowManagerSystem {
  protected def systemName = "cromwell-system"

  protected def newActorSystem(): ActorSystem = ActorSystem(systemName)

  implicit final lazy val actorSystem = newActorSystem()

  def shutdownActorSystem(): Unit = {
    actorSystem.shutdown()
  }

  CromwellBackend.initBackends(BackendConfiguration.AllBackendEntries, BackendConfiguration.DefaultBackendEntry, actorSystem)
  // For now there's only one WorkflowManagerActor so no need to dynamically name it
  lazy val workflowManagerActor = actorSystem.actorOf(WorkflowManagerActor.props(), "WorkflowManagerActor")
}