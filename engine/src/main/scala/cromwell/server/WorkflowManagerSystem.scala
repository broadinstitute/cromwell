
package cromwell.server

import akka.actor.ActorSystem
import com.typesafe.config.ConfigFactory
import cromwell.engine.backend.{BackendConfiguration, CromwellBackend}
import cromwell.engine.workflow.{ShadowWorkflowManagerActor, WorkflowManagerActor}

trait WorkflowManagerSystem {
  protected def systemName = "cromwell-system"

  protected def newActorSystem(): ActorSystem = ActorSystem(systemName)
  val conf = ConfigFactory.load()
  val shadowMode = conf.getBoolean("system.shadowExecutionEnabled")

  implicit final lazy val actorSystem = newActorSystem()

  def shutdownActorSystem(): Unit = {
    actorSystem.shutdown()
  }

  CromwellBackend.initBackends(
    BackendConfiguration.AllBackendEntries,
    BackendConfiguration.DefaultBackendEntry,
    actorSystem,
    ConfigFactory.load.getBoolean("system.shadowExecutionEnabled")
  )

  // For now there's only one WorkflowManagerActor so no need to dynamically name it
  val workflowManagerProps = if(shadowMode) ShadowWorkflowManagerActor.props() else WorkflowManagerActor.props()
  lazy val workflowManagerActor = actorSystem.actorOf(workflowManagerProps, "WorkflowManagerActor")
}