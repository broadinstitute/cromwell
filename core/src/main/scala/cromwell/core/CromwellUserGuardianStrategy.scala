package cromwell.core

import akka.actor.SupervisorStrategy.Escalate
import akka.actor.{ActorInitializationException, OneForOneStrategy, SupervisorStrategy, SupervisorStrategyConfigurator}

class CromwellUserGuardianStrategy extends SupervisorStrategyConfigurator {
  override def create(): SupervisorStrategy = OneForOneStrategy() {
    case _: ActorInitializationException => Escalate
    case t => akka.actor.SupervisorStrategy.defaultDecider.applyOrElse(t, (_: Any) => Escalate)
  }
}
