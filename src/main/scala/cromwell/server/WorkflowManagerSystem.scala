package cromwell.server

import akka.actor.ActorSystem
import com.typesafe.config.ConfigFactory
import cromwell.engine.backend.{Backend, CromwellBackend}
import cromwell.engine.workflow.WorkflowManagerActor

object GlobalActorSystem {
  val globalSystemName = "cromwell-system"
  final lazy val globalActorSystem = ActorSystem(globalSystemName)
}

trait WorkflowManagerSystem {
  import GlobalActorSystem._

  protected def systemName = globalSystemName
  protected def newActorSystem(): ActorSystem = globalActorSystem
  implicit val actorSystem = newActorSystem()

  def shutdownActorSystem(): Unit = {
    actorSystem.shutdown()
  }

  def backendType: String = ConfigFactory.load.getConfig("backend").getString("backend")

  lazy val backend: Backend = CromwellBackend.initBackend(backendType, actorSystem)
  // For now there's only one WorkflowManagerActor so no need to dynamically name it
  lazy val workflowManagerActor = actorSystem.actorOf(WorkflowManagerActor.props(backend), "WorkflowManagerActor")
}
