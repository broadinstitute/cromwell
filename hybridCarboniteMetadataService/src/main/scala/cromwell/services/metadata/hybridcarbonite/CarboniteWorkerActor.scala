package cromwell.services.metadata.hybridcarbonite

import akka.actor.{Actor, ActorLogging, ActorRef, Props}
import cromwell.util.GracefulShutdownHelper
import cromwell.util.GracefulShutdownHelper.ShutdownCommand
import cromwell.services.metadata.hybridcarbonite.CarboniteWorkerActor.DoCarboniting

import scala.concurrent.ExecutionContext
import scala.concurrent.duration._

class CarboniteWorkerActor(carboniteInterval: FiniteDuration, ioActor: ActorRef) extends Actor with ActorLogging with GracefulShutdownHelper {

  implicit val ec: ExecutionContext = context.dispatcher

  scheduleNextCarboniting()

  override def receive: Receive = {
    case DoCarboniting => doCarboniting()
    case ShutdownCommand => stopGracefully()
    case oops => log.error(s"Programmer Error! The CarboniteWorkerActor is not to be talked to! ($sender sent $oops})")
  }

  def doCarboniting(): Unit = {
    // TODO: [CARBONITE] When carboniting actually does something, we probably metrics here rather than a tick log...
    log.info("Carbonite Worker Tick...")

    scheduleNextCarboniting()
  }

  def scheduleNextCarboniting(): Unit = {
    context.system.scheduler.scheduleOnce(carboniteInterval) {
      self ! DoCarboniting
    }
    ()
  }

  // TODO: [CARBONITE] When the carboniting process is implemented, we might need to implement graceful shutdowns
  def stopGracefully(): Unit = {
    context.stop(self)
  }
}

object CarboniteWorkerActor {
  def props(carboniteInterval: FiniteDuration, ioActor: ActorRef) = Props(new CarboniteWorkerActor(carboniteInterval, ioActor))

  case object DoCarboniting
}
