package cromwell.engine

import akka.actor.{Actor, Props}
import akka.event.{Logging, LoggingReceive}
import cromwell.engine.backend._

import scala.language.postfixOps

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

  private def execute(): Unit = {
    log.info(s"$tag: starting ${backendCall.call.name} for workflow ${backendCall.workflowDescriptor.shortId}")
    val result = backendCall.execute

    result match {
      case SuccessfulExecution(_) => log.info(s"$tag: successful execution.")
      case AbortedExecution => log.info(s"$tag: aborted.")
      case FailedExecution(e, returnCode) => log.error(e, s"$tag: failed execution, returnCode = $returnCode")
    }

    context.parent ! CallActor.ExecutionFinished(backendCall.call, result)
    context.stop(self)
  }

  override def receive = LoggingReceive {
    case Execute => execute()
    case badMessage => log.error(s"$tag: unexpected message $badMessage.")
  }
}
