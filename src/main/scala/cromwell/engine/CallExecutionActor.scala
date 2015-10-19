package cromwell.engine

import akka.actor.{Actor, Props}
import akka.event.{Logging, LoggingReceive}
import cromwell.engine.backend._

import scala.language.postfixOps
import scala.util.Try

object CallExecutionActor {
  sealed trait CallExecutionActorMessage

  sealed trait ExecutionMode extends CallExecutionActorMessage {
    def execute(backendCall: BackendCall): ExecutionResult
  }

  case object Execute extends ExecutionMode {
    override def execute(backendCall: BackendCall) = backendCall.execute
  }

  final case class Resume(jobKey: JobKey) extends ExecutionMode {
    override def execute(backendCall: BackendCall) = backendCall.resume(jobKey)
  }

  def props(backendCall: BackendCall): Props = Props(new CallExecutionActor(backendCall))
}

/** Actor to manage the execution of a single call. */
class CallExecutionActor(backendCall: BackendCall) extends Actor with CromwellActor {
  import CallExecutionActor._

  private val log = Logging(context.system, classOf[CallExecutionActor])

  val tag = s"CallExecutionActor [UUID(${backendCall.workflowDescriptor.shortId}):${backendCall.key.tag}]"

  private def execute(executionMode: ExecutionMode): Unit = {
    log.info(s"$tag: starting ${backendCall.call.name} for workflow ${backendCall.workflowDescriptor.shortId}")
    // If the actual execution throws an exception, catch and wrap with a `FailedExecution`.
    // Ideally the backend would map a command execution failure to a `FailedExecution`, but
    // if an exception propagates out of `execute` it should not take down this actor.
    val result = Try(executionMode.execute(backendCall)) recover { case e: Exception => FailedExecution(e, None) } get

    result match {
      case SuccessfulExecution(_, _) => log.info(s"$tag: successful execution.")
      case AbortedExecution => log.info(s"$tag: aborted.")
      case FailedExecution(e, returnCode) => log.error(e, s"$tag: failed execution, returnCode = $returnCode")
    }

    context.parent ! CallActor.ExecutionFinished(backendCall.call, result)
    context.stop(self)
  }

  override def receive = LoggingReceive {
    case mode: ExecutionMode => execute(mode)
    case badMessage => log.error(s"$tag: unexpected message $badMessage.")
  }
}
