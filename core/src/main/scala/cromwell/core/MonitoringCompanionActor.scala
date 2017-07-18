package cromwell.core

import akka.actor.{Actor, ActorLogging, ActorRef, Props}
import cromwell.core.MonitoringCompanionActor._
import cromwell.util.GracefulShutdownHelper.ShutdownCommand

import scala.concurrent.duration._
import scala.language.postfixOps

object MonitoringCompanionActor {
  sealed trait MonitoringCompanionCommand
  private [core] case object AddWork extends MonitoringCompanionCommand
  private [core] case object RemoveWork extends MonitoringCompanionCommand
  private [core] def props(actorToMonitor: ActorRef) = Props(new MonitoringCompanionActor(actorToMonitor))
}

private [core] class MonitoringCompanionActor(actorToMonitor: ActorRef) extends Actor with ActorLogging {
  private var workCount: Int = 0
  
  override def receive = {
    case AddWork => workCount += 1
    case RemoveWork => workCount -= 1
    case ShutdownCommand if workCount <= 0 =>
      context stop actorToMonitor
      context stop self
    case ShutdownCommand => 
      log.info(s"{} is still processing {} messages", actorToMonitor.path.name, workCount)
      context.system.scheduler.scheduleOnce(1 second, self, ShutdownCommand)(context.dispatcher)
      ()
  }
}

trait MonitoringCompanionHelper { this: Actor =>
  private val monitoringActor = context.actorOf(MonitoringCompanionActor.props(self))
  private var shuttingDown: Boolean = false
  
  def addWork() = monitoringActor ! AddWork
  def removeWork() = monitoringActor ! RemoveWork

  val monitoringReceive: Receive = {
    case ShutdownCommand if !shuttingDown => 
      shuttingDown = true
      monitoringActor ! ShutdownCommand
    case ShutdownCommand => // Ignore if we're already shutting down
  }
}
