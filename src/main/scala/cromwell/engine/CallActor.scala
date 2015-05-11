package cromwell.engine

import akka.actor.{ActorRef, Actor}
import cromwell.binding.{Call, WdlValue}
import cromwell.engine.CallActor.{Done, Started, Run, Start}
import scala.language.postfixOps
import sys.process._

object CallActor {

  case class Start(call: Call, parameters: Map[String, WdlValue])

  case class Run(commandLine: String, workflowActor: ActorRef)

  case object Started

  case object Done

}

/**
 * Probably not the right level of abstraction for a system that supports richer
 * control structures like loops, conditionals, scatter etc., but this should be
 * sufficient for Sprint 2.
 */
class CallActor extends Actor {

  override def receive: Receive = {

    case Start(call, parameters) =>
      // Assuming there is no check required here for presence and type of
      // parameters, that can be added if needed.
      val commandLine = call.task.command.instantiate(parameters)
      self ! Run(commandLine, sender())
      sender ! Started

    case Run(commandLine, workflowActor) =>
      // TODO This is an "EchoBackend", need to decouple from a particular backend
      s"echo $commandLine" !!;
      workflowActor ! Done
  }
}
