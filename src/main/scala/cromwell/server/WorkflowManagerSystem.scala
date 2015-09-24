package cromwell.server

import akka.actor.ActorSystem
import cromwell.engine.db.DataAccess
import cromwell.engine.workflow.WorkflowManagerActor
import cromwell.parser.BackendType

trait WorkflowManagerSystem {
  lazy val backendType: BackendType = WorkflowManagerActor.BackendType

  val systemName = "cromwell-system"
  implicit val actorSystem = ActorSystem(systemName)

  // For now there's only one WorkflowManagerActor so no need to dynamically name it
  lazy val workflowManagerActor = actorSystem.actorOf(WorkflowManagerActor.props(dataAccess, backendType), "WorkflowManagerActor")

  // Lazily created as the primary consumer is the workflowManagerActor.
  protected lazy val dataAccess: DataAccess = DataAccess()

  /**
   * Should be called after the system is no longer in use.
   *
   * @return Non-blocking future that will eventually shut down this instance or return an error.
   */
  def shutdown() = dataAccess.shutdown()
}
