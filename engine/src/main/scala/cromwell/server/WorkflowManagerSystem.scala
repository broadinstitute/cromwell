package cromwell.server

import akka.actor.ActorSystem
import com.typesafe.config.ConfigFactory
import cromwell.engine.backend.CromwellBackend
import cromwell.engine.workflow.WorkflowManagerActor
import collection.JavaConversions._

trait WorkflowManagerSystem {
  protected def systemName = "cromwell-system"

  protected def newActorSystem(): ActorSystem = ActorSystem(systemName)
  val conf = ConfigFactory.load()
  val shadowMode = conf.getBoolean("system.shadowExecutionEnabled")

  implicit final lazy val actorSystem = newActorSystem()

  def shutdownActorSystem(): Unit = {
    actorSystem.shutdown()
  }

  def allowedBackends: List[String] = ConfigFactory.load.getConfig("backend").getStringList("backendsAllowed").toList

  def defaultBackend: String = ConfigFactory.load.getConfig("backend").getString("defaultBackend")

  if (!shadowMode) CromwellBackend.initBackends(allowedBackends, defaultBackend, actorSystem)
  // For now there's only one WorkflowManagerActor so no need to dynamically name it
  lazy val workflowManagerActor = actorSystem.actorOf(WorkflowManagerActor.props(shadowMode), "WorkflowManagerActor")
}
