package cromwell.engine

import akka.actor.{Actor, Props}
import akka.event.{Logging, LoggingReceive}
import cromwell.binding.WdlExpression.ScopedLookupFunction
import cromwell.binding.{Call, CallInputs, WorkflowDescriptor}
import cromwell.engine.backend.{BackendCall, TaskAbortedException, Backend}

import scala.language.postfixOps
import scala.util.{Success, Failure}

object CallExecutionActor {
  sealed trait CallExecutionActorMessage
  case object Execute extends CallExecutionActorMessage
  def props(backendCall: BackendCall): Props = Props(new CallExecutionActor(backendCall))
}

/** Actor to manage the execution of a single call. */
class CallExecutionActor(backendCall: BackendCall) extends Actor with CromwellActor {
  import CallExecutionActor._

  private val log = Logging(context.system, classOf[CallExecutionActor])

  def tag = s"CallExecutionActor [UUID(${backendCall.workflowDescriptor.shortId}):${backendCall.call.name}]"

  private def execute = {
    log.info(s"$tag: starting ${backendCall.call.name} for workflow ${backendCall.workflowDescriptor.shortId}")
    val result = backendCall.execute
    result match {
      case Success(_) => log.info(s"$tag: successful execution.")
      case Failure(e: TaskAbortedException) => log.info(s"$tag: aborted.")
      case Failure(e) => log.error(s"$tag: failed execution - ${e.getMessage}")
    }
    context.parent ! CallActor.ExecutionFinished(backendCall.call, result)
  }

  override def receive = LoggingReceive {
    case Execute => execute
    case badMessage => log.error(s"$tag: unexpected message $badMessage.")
  }
}
