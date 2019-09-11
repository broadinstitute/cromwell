package cromwell.services.metadata.hybridcarbonite

import akka.actor.{Actor, ActorLogging, Props}
import cromwell.util.GracefulShutdownHelper
import cromwell.util.GracefulShutdownHelper.ShutdownCommand
import cromwell.services.metadata.hybridcarbonite.CarboniteWorkerActor.DoCarboniting

import scala.concurrent.ExecutionContext
import scala.concurrent.duration._

class CarboniteWorkerActor(carboniteInterval: FiniteDuration) extends Actor with ActorLogging with GracefulShutdownHelper {

  implicit val ec: ExecutionContext = context.dispatcher

  scheduleNextCarboniting()

  override def receive: Receive = {
    case DoCarboniting => doCarboniting()
    case ShutdownCommand => stopGracefully()
    case oops => log.error(s"Programmer Error! The CarboniteWorkerActor is not to be talked to! ($sender sent $oops})")
  }

  def doCarboniting(): Unit = {
    log.info("Carbonite Worker Tick...")
    scheduleNextCarboniting()
  }

  def scheduleNextCarboniting(): Unit = {
    context.system.scheduler.scheduleOnce(carboniteInterval) {
      self ! DoCarboniting
    }
    ()
  }

  // TODO: Can we be more graceful?
  def stopGracefully(): Unit = {
    context.stop(self)
  }
}

object CarboniteWorkerActor {
  def props(carboniteInterval: FiniteDuration) = Props(new CarboniteWorkerActor(carboniteInterval))

  case object DoCarboniting
}
