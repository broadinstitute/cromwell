package cromwell.util

import akka.actor.SupervisorStrategy.{Decider, Stop}
import akka.actor.{Actor, ActorRef, OneForOneStrategy, SupervisorStrategy}

trait StopAndLogSupervisor { this: Actor =>

  private var failureLog: Map[ActorRef, Throwable] = Map.empty

  final val stopAndLogStrategy: SupervisorStrategy = {
    def stoppingDecider: Decider = {
      case e: Exception =>
        val failer = sender()
        failureLog += failer -> e
        Stop
    }
    OneForOneStrategy(loggingEnabled = false)(stoppingDecider)
  }

  final def getFailureCause(actorRef: ActorRef): Option[Throwable] = {
    val result = failureLog.get(actorRef)
    failureLog -= actorRef
    result
  }

  override final val supervisorStrategy = stopAndLogStrategy
}
