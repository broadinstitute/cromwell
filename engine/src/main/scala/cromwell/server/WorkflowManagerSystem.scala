
package cromwell.server

import akka.actor.ActorSystem
import com.typesafe.config.ConfigFactory
import cromwell.engine.backend.{BackendConfiguration, CromwellBackends}
import cromwell.engine.workflow.workflowstore.{FileSystemWorkflowStore$, WorkflowStoreActor}
import cromwell.engine.workflow.WorkflowManagerActor
import org.slf4j.LoggerFactory

trait WorkflowManagerSystem {
  protected def systemName = "cromwell-system"

  protected def newActorSystem(): ActorSystem = ActorSystem(systemName)
  val conf = ConfigFactory.load()
  val logger = LoggerFactory.getLogger(getClass.getName)
  implicit final lazy val actorSystem = newActorSystem()

  def shutdownActorSystem(): Unit = {
    actorSystem.shutdown()
  }

  CromwellBackends.initBackends(
    BackendConfiguration.AllBackendEntries,
    actorSystem
  )

  lazy val workflowStoreActor = actorSystem.actorOf(WorkflowStoreActor.props(FileSystemWorkflowStore$), "WorkflowStoreActor")

  // For now there's only one WorkflowManagerActor so no need to dynamically name it
  val workflowManagerProps = WorkflowManagerActor.props(workflowStoreActor)
  lazy val workflowManagerActor = actorSystem.actorOf(workflowManagerProps, "WorkflowManagerActor")
}
