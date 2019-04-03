package cromwell.util

import akka.actor.{Actor, ActorRef, ActorSystem, Kill, PoisonPill, SupervisorStrategy}

import scala.util.control.NoStackTrace

object AkkaTestUtil {

  def actorDeathMethods(system: ActorSystem): List[(String, ActorRef => Unit)] = List(
    ("external_stop", (a: ActorRef) => system.stop(a)),
    ("internal_stop", (a: ActorRef) => a ! InternalStop),
    ("poison_pill", (a: ActorRef) => a ! PoisonPill),
    ("kill_message", (a: ActorRef) => a ! Kill),
    ("throw_exception", (a: ActorRef) => a ! ThrowException)
  )

  case object InternalStop
  case object ThrowException

  class StoppingSupervisor extends Actor {
    override val supervisorStrategy = SupervisorStrategy.stoppingStrategy
    def receive = Actor.emptyBehavior
  }

  class DeathTestActor extends Actor {
    protected def stoppingReceive: Actor.Receive = {
      case InternalStop => context.stop(self)
      case ThrowException => throw new Exception("Don't panic, dear debugger! This was a deliberate exception for the test case.") with NoStackTrace
    }
    override def receive = stoppingReceive orElse Actor.ignoringBehavior
  }
}
