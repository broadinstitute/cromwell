package cromwell.engine

import akka.actor.{Actor, ActorRef, Props}
import akka.event.{Logging, LoggingReceive}
import cromwell.binding.Call
import cromwell.binding.values.WdlValue
import cromwell.engine.CallActor._
import cromwell.engine.WorkflowActor.InvalidOperation
import cromwell.engine.backend.Backend
import cromwell.util.TryUtil

import scala.language.postfixOps
import scala.util.Success


object CallActor {

  sealed trait CallActorMessage

  case class Start(call: Call, backend: Backend, symbolStore: SymbolStore) extends CallActorMessage
  // Does this message need to pass along the SymbolStore for output expression evaluation?
  case class Run(commandLine: String, backend: Backend, call: Call, symbolStore: SymbolStore) extends CallActorMessage
  case class Started(call: Call) extends CallActorMessage
  case class Done(call: Call, outputs: Map[String, WdlValue]) extends CallActorMessage
  case class Failed(call: Call, failure: String) extends CallActorMessage

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

  private def notStarted: Receive = {
    LoggingReceive {
      case Start(call, backend, symbolStore) =>
        val inputs = symbolStore.locallyQualifiedInputs(call)
        val tryCommandLine = call.task.command.instantiate(inputs)
        tryCommandLine match {
          case Success(commandLine) =>
            _creator = Option(sender())
            self ! Run(commandLine, backend, call, symbolStore)
            creator ! Started(call)
            context.become(started)
          case failure =>
            sender ! failure
        }

      case bad@_ if !isStarted =>
        sender ! InvalidOperation(s"Received message $bad, but not yet Started.")
    }
  }

  private def started: Receive = {
    LoggingReceive {
      case Run(commandLine, backend, call, symbolStore) =>
        log.info(self + " executing command: " + commandLine)
        val tryOutputs = backend.executeCommand(commandLine, call, call.task.outputs, symbolStore)
        val (successes, failures) = tryOutputs.partition { _._2.isSuccess }

        if (failures.isEmpty) {
          // Materialize the Successes.
          val outputs = successes.map { case (key, value) => key -> value.get }
          creator ! Done(call, outputs)
        } else {
          val errorMessages = TryUtil.stringifyFailures(failures.values).mkString("\n")

          log.error(errorMessages)
          creator ! Failed(call, errorMessages)
        }

      case unknown @_ =>
        sender ! InvalidOperation(s"Unexpected message $unknown.")
    }
  }

  override def receive = notStarted
}
