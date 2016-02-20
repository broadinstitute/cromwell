package cromwell.server

import akka.actor.ActorSystem
import com.typesafe.config.ConfigFactory
import cromwell.engine.workflow.WorkflowManagerActor

// TODO: WARNING! This absolutely *NEEDS* to be gone. It's here just so that it can complement CromwellBackend singleton object, which in turn *must* also be gone!
object WorkflowManagerSystem {
  //This is just like an accessor to get the ActorSystem
  var system: ActorSystem = _
  def setActorSystem(actorSystem: ActorSystem): Unit = {
    system = actorSystem
  }
}

trait WorkflowManagerSystem {
  protected def systemName = "cromwell-system"
  protected def newActorSystem(): ActorSystem = ActorSystem(systemName)
  implicit final val actorSystem = newActorSystem()
  WorkflowManagerSystem.setActorSystem(actorSystem)

  def shutdownActorSystem(): Unit = {
    actorSystem.shutdown()
  }

  //  lazy val backend: Backend = CromwellBackend.initBackend(backendType, actorSystem)
  // For now there's only one WorkflowManagerActor so no need to dynamically name it
  lazy val workflowManagerActor = actorSystem.actorOf(WorkflowManagerActor.props(), "WorkflowManagerActor")
}
