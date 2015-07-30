package cromwell.engine

import akka.actor.{Actor, Props}
import akka.event.{Logging, LoggingReceive}
import cromwell.binding.WdlExpression.ScopedLookupFunction
import cromwell.binding.{Call, CallInputs, WorkflowDescriptor}
import cromwell.engine.backend.{TaskAbortedException, Backend}

import scala.language.postfixOps
import scala.util.{Success, Failure}

object CallExecutionActor {

  sealed trait CallExecutionActorMessage
  case class Execute(workflowId: WorkflowId,
                     backend: Backend,
                     command: String,
                     workflowDescriptor: WorkflowDescriptor,
                     call: Call,
                     callInputs: CallInputs,
                     lookup: ScopedLookupFunction) extends CallExecutionActorMessage

  def props(callReference: CallReference): Props = Props(new CallExecutionActor(callReference))
}

/** Actor to manage the execution of a single call. */
class CallExecutionActor(callReference: CallReference) extends Actor with CromwellActor {
  import CallExecutionActor._

  private val log = Logging(context.system, classOf[CallExecutionActor])

  def tag = s"CallExecutionActor [$callReference]"

  private def execute(workflowId: WorkflowId,
                        backend: Backend,
                        command: String,
                        workflowDescriptor: WorkflowDescriptor,
                        call: Call,
                        callInputs: CallInputs,
                        lookup: ScopedLookupFunction) = {

    log.info(s"$tag: starting ${call.name} for workflow ${workflowDescriptor.shortId}")

    val executionResult = backend.executeCommand(
      command,
      workflowDescriptor,
      call,
      callInputs,
      lookup,
      AbortFunctionRegistration(registerAbortFunction(callReference)))

    executionResult match {
      case Success(_) => log.info(s"$tag: successful execution.")
      case Failure(e: TaskAbortedException) => log.info(s"$tag: aborted.")
      case Failure(e) => log.error(s"$tag: failed execution - ${e.getMessage}")
    }

    context.parent ! CallActor.ExecutionFinished(call, executionResult)
  }

  override def receive = LoggingReceive {
    case Execute(workflowId, backend, command, workflowDescription, call, callInputs, lookup) =>
      execute(workflowId, backend, command, workflowDescription, call, callInputs, lookup)
    case badMessage =>
      val diagnostic = s"$tag: unexpected message $badMessage."
      log.error(diagnostic)
  }
  
  private def registerAbortFunction(callReference: CallReference)(abortFunction: AbortFunction): Unit = {
    context.parent ! CallActor.RegisterCallAbortFunction(abortFunction)
  }
}
