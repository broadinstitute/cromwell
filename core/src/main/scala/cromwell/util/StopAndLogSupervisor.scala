package cromwell.util

import akka.actor.SupervisorStrategy.{Decider, Stop}
import akka.actor.{Actor, ActorRef, OneForOneStrategy, SupervisorStrategy}
import cromwell.core.logging.WorkflowLogging

trait StopAndLogSupervisor { this: Actor with WorkflowLogging =>

  private var failureLog: Map[ActorRef, Throwable] = Map.empty

  final val stopAndLogStrategy: SupervisorStrategy = {
    def stoppingDecider: Decider = {
      case e: Exception =>
        val failer = sender()
        failureLog += failer -> e
        Stop
    }
    OneForOneStrategy()(stoppingDecider)
  }

  final def getFailureCause(actorRef: ActorRef): Option[Throwable] = failureLog.get(actorRef)

  override final val supervisorStrategy = stopAndLogStrategy
}
