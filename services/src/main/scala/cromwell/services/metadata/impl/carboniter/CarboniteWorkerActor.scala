package cromwell.services.metadata.impl.carboniter

import akka.actor.{Actor, ActorLogging, Props}
import cromwell.services.metadata.impl.carboniter.CarboniteWorkerActor.DoCarboniting
import scala.concurrent.ExecutionContext
import scala.concurrent.duration._

class CarboniteWorkerActor(carboniteInterval: FiniteDuration) extends Actor with ActorLogging {

  implicit val ec: ExecutionContext = context.dispatcher

  scheduleNextCarboniting()

  override def receive: Receive = {
    case DoCarboniting => doCarboniting()

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
}

object CarboniteWorkerActor {
  def props(carboniteInterval: FiniteDuration) = Props(new CarboniteWorkerActor(carboniteInterval))

  case object DoCarboniting
}
