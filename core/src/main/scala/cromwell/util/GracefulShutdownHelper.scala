package cromwell.util

import akka.actor.{Actor, ActorLogging, ActorRef, Terminated}
import akka.pattern.GracefulStopSupport
import cats.data.NonEmptyList
import cromwell.util.GracefulShutdownHelper.ShutdownCommand

object GracefulShutdownHelper {
  case object ShutdownCommand
}

trait GracefulShutdownHelper extends GracefulStopSupport { this: Actor with ActorLogging =>
  private var shuttingDown: Boolean = false
  private var shutdownList: Set[ActorRef] = Set.empty
  
  def isShuttingDown: Boolean = shuttingDown

  def waitForActorsAndShutdown(actorsLists: NonEmptyList[ActorRef]): Unit = {
    if (shuttingDown) {
      log.error("Programmer error, this actor has already initiated its shutdown. Only call this once per actor !")
    } else {
      shuttingDown = true
      shutdownList = actorsLists.toList.toSet
      shutdownList foreach context.watch
      shutdownList foreach { _ ! ShutdownCommand }
      
      context become {
        case Terminated(actor) if shuttingDown && shutdownList.contains(actor) =>
          shutdownList = shutdownList - actor
          if (shutdownList.isEmpty) context stop self
      }
    }
  }
}
