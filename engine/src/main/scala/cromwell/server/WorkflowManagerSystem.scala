
package cromwell.server

import akka.actor.ActorSystem
import com.typesafe.config.ConfigFactory
import cromwell.engine.backend.{BackendConfiguration, CromwellBackend}
import cromwell.engine.workflow.{ShadowWorkflowManagerActor, WorkflowManagerActor}
import cromwell.logging.WorkflowLogger
import org.slf4j.LoggerFactory

trait WorkflowManagerSystem {
  protected def systemName = "cromwell-system"

  protected def newActorSystem(): ActorSystem = ActorSystem(systemName)
  val conf = ConfigFactory.load()
  val logger = LoggerFactory.getLogger(getClass.getName)

  val mulletModeOption = "system.mulletMode"
  lazy val mulletMode = if (conf.hasPath(mulletModeOption)) conf.getBoolean(mulletModeOption) else false

  // TODO: At some point we won't have a retro mode anymore. In the meantime, you have to actively choose to be retro:
  if (conf.hasPath("system.shadowExecutionEnabled")) logger.warn("Shadow mode is now the default. You must now actively choose to be retro by setting system.mulletMode = true")

  implicit final lazy val actorSystem = newActorSystem()

  def shutdownActorSystem(): Unit = {
    actorSystem.shutdown()
  }

  if (mulletMode) CromwellBackend.initBackends(BackendConfiguration.AllBackendEntries, BackendConfiguration.DefaultBackendEntry, actorSystem)
  // For now there's only one WorkflowManagerActor so no need to dynamically name it
  val workflowManagerProps = if (mulletMode) WorkflowManagerActor.props() else ShadowWorkflowManagerActor.props()
  lazy val workflowManagerActor = actorSystem.actorOf(workflowManagerProps, "WorkflowManagerActor")
}