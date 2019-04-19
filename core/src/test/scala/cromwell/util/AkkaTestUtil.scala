package cromwell.util

import akka.actor.{Actor, ActorLogging, ActorRef, ActorSystem, Kill, PoisonPill, Props, SupervisorStrategy}
import akka.testkit.TestProbe

import scala.util.control.NoStackTrace

object AkkaTestUtil {

  implicit class EnhancedTestProbe(probe: TestProbe) {

    // Get a 'props' handle which will return a (inbound-only!!) wrapper around this test probe when constructed
    def props = Props(new Actor with ActorLogging {
      def receive = {
        case outbound @ _ if sender == probe.ref =>
          val msg = "Unexpected outbound message from Probe. You're doing something wrong!"
          log.error(msg)
          throw new RuntimeException(msg)
        case inbound => probe.ref forward inbound
      }
    })
  }

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
