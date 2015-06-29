package cromwell.server

import akka.actor.{ActorRef, ActorSystem}
import com.typesafe.config.ConfigFactory
import cromwell.engine.WorkflowManagerActor

trait WorkflowManagerSystem {
  val systemName: String
  implicit val actorSystem: ActorSystem
  val workflowManagerActor: ActorRef
}
  
case class DefaultWorkflowManagerSystem() extends WorkflowManagerSystem {
  val systemName = "cromwell-system"
  implicit val actorSystem = ActorSystem(systemName)
  actorSystem.registerOnTermination {actorSystem.log.info(s"$systemName shutting down")}
  val workflowManagerActor = actorSystem.actorOf(WorkflowManagerActor.props)
}
