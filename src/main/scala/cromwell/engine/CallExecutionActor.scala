package cromwell.engine

import akka.actor.{Actor, Props}
import akka.event.{Logging, LoggingReceive}
import cromwell.binding.WdlExpression.ScopedLookupFunction
import cromwell.binding.{Call, CallInputs, WorkflowDescriptor}
import cromwell.engine.backend.Backend

import scala.language.postfixOps
import scala.util.{Failure, Success}


object CallExecutionActor {
  def props(backend: Backend, command: String, workflowDescriptor: WorkflowDescriptor, call: Call, callInputs: CallInputs, lookup: ScopedLookupFunction): Props =
    Props(new CallExecutionActor(backend, command, workflowDescriptor, call, callInputs, lookup))
}

/** Actor to manage the execution of a single call. */
class CallExecutionActor(backend: Backend, command: String, workflowDescriptor: WorkflowDescriptor, call: Call, callInputs: CallInputs, lookup: ScopedLookupFunction) extends Actor with CromwellActor {
  private val log = Logging(context.system, classOf[CallExecutionActor])
  val tag = s"CallExecutionActor [UUID(${workflowDescriptor.shortId}):${call.name}]"

  log.info(s"$tag: starting.")
  backend.executeCommand(command, workflowDescriptor, call, callInputs, lookup) match {
    case Success(callOutputs) =>
      log.info(s"$tag: successful execution.")
      context.parent ! CallActor.ExecutionFinished(callOutputs)
    case Failure(e) =>
      log.error(s"$tag: failed execution.")
      context.parent ! CallActor.ExecutionFailed(e)
  }

  /**
   * The `CallExecutionActor` is not meant to respond to messages, it just starts up, does its work,
   * and shuts down.
   */
  override def receive = LoggingReceive {
    case badMessage =>
      val diagnostic = s"$tag: unexpected message $badMessage."
      log.error(diagnostic)
  }
}
