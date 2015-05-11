package cromwell.engine

import akka.actor.SupervisorStrategy.Stop
import akka.actor.{Actor, ActorRef, Props}
import akka.event.{Logging, LoggingReceive}
import cromwell.binding.{Call, WdlValue}
import cromwell.engine.CallActor._
import cromwell.engine.WorkflowActor.InvalidOperation

import scala.language.postfixOps
import scala.util.{Failure, Success}


object CallActor {

  sealed trait CallActorMessage

  case class Start(call: Call, symbolStore: SymbolStore) extends CallActorMessage
  case class Run(commandLine: String) extends CallActorMessage
  case object Started extends CallActorMessage
  case class Done(outputs: Map[String, WdlValue]) extends CallActorMessage
  case class Failed(throwable: Throwable) extends CallActorMessage
  case object Stopped extends CallActorMessage

  def props: Props = Props(classOf[CallActor])
}

/**
 * Probably not the right level of abstraction for a system that supports richer
 * control structures like loops, conditionals, scatter etc., but this should be
 * sufficient for Sprint 2/3.
 */
class CallActor extends Actor {

  private val log = Logging(context.system, this)

  private var _creator: Option[ActorRef] = None

  private def isStarted: Boolean = _creator.isDefined

  private def creator: ActorRef = _creator.get

  override def receive: Receive =
    LoggingReceive {
      case Start(call, symbolStore) =>
        val inputs = symbolStore.locallyQualifiedInputs(call)
        val tryCommandLine = call.task.command.instantiate(inputs)
        tryCommandLine match {
          case Success(commandLine) =>
            _creator = Option(sender())
            self ! Run(commandLine)
            creator ! Started
          case Failure(t) =>
            sender ! Failed(t)
        }

      case bad @_ if !isStarted =>
        sender ! InvalidOperation(s"Received message $bad, but not yet Started.")

      case Run(commandLine) =>
        log.info(commandLine)
        // TODO get actual outputs and send them along in the reply
        creator ! Done(Map.empty[String, WdlValue])
        self ! Stop

      case Stop =>
        context.stop(self)
        creator ! Stopped

      case unknown @_ =>
        sender ! InvalidOperation(s"Unexpected message $unknown.")
    }
}
