package cromwell.server

import akka.actor.{Props, ActorSystem}
import com.typesafe.config.ConfigFactory
import cromwell.engine.backend.CromwellBackend
import cromwell.engine.workflow.{MaterializeWorkflowDescriptorActor, WorkflowManagerActor}
import cromwell.engine.workflow.WorkflowManagerActor
import cromwell.instrumentation.Instrumentation._
import collection.JavaConversions._

trait WorkflowManagerSystem {
  Monitor.start()

  protected def systemName = "cromwell-system"

  protected def newActorSystem(): ActorSystem = ActorSystem(systemName)

  implicit final lazy val actorSystem = newActorSystem()

  def shutdownActorSystem(): Unit = {
    actorSystem.shutdown()
  }

  def allowedBackends: List[String] = ConfigFactory.load.getConfig("backend").getStringList("backendsAllowed").toList
  def defaultBackend: String = ConfigFactory.load.getConfig("backend").getString("defaultBackend")

  CromwellBackend.initBackends(allowedBackends, defaultBackend, actorSystem)
  // For now there's only one WorkflowManagerActor so no need to dynamically name it
  lazy val workflowManagerActor = actorSystem.actorOf(WorkflowManagerActor.props(), "WorkflowManagerActor")
}
