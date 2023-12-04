package cromwell.util

import akka.actor.SupervisorStrategy.{Decider, Stop}
import akka.actor.{Actor, ActorRef, OneForOneStrategy, SupervisorStrategy}

trait StopAndLogSupervisor { this: Actor =>

  protected def onFailure(actorRef: ActorRef, throwable: => Throwable): Unit

  final val stopAndLogStrategy: SupervisorStrategy = {
    def stoppingDecider: Decider = { case e: Exception =>
      onFailure(sender(), e)
      Stop
    }
    OneForOneStrategy(loggingEnabled = false)(stoppingDecider)
  }

  final override val supervisorStrategy = stopAndLogStrategy
}
