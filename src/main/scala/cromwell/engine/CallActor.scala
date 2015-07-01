package cromwell.engine

import akka.actor.{Actor, Props}
import akka.event.{Logging, LoggingReceive}
import cromwell.binding.values.WdlValue
import cromwell.binding.{Call, CallInputs, WorkflowDescriptor}
import cromwell.engine.backend.Backend
import cromwell.engine.workflow.WorkflowActor
import cromwell.engine.workflow.WorkflowActor.CallFailed
import cromwell.util.TryUtil

import scala.language.postfixOps
import scala.util.{Failure, Success}


object CallActor {
  sealed trait CallActorMessage
  case class Start(inputs: CallInputs) extends CallActorMessage

  def props(call: Call, locallyQualifiedInputs: Map[String, WdlValue], backend: Backend, workflowDescriptor: WorkflowDescriptor): Props =
    Props(new CallActor(call, locallyQualifiedInputs, backend, workflowDescriptor))
}


/** Actor to manage the execution of a single call. */
class CallActor(call: Call, locallyQualifiedInputs: Map[String, WdlValue], backend: Backend, workflowDescriptor: WorkflowDescriptor) extends Actor with CromwellActor {

  private val log = Logging(context.system, classOf[CallActor])
  val tag = s"CallActor [UUID(${workflowDescriptor.shortId}):${call.name}]"

  override def receive = LoggingReceive {
    case CallActor.Start => handleStart()
    case badMessage =>
      val diagnostic = s"$tag: unexpected message $badMessage."
      log.error(diagnostic)
      context.parent ! WorkflowActor.CallFailed(call, diagnostic)
  }

  /**
   * Performs the following steps:
   * <ol>
   *   <li>Instantiates the command line using these inputs.</li>
   *   <li>If the above completes successfully, messages the sender with `Started`,
   *   otherwise messages `Failed`.</li>
   *   <li>Executes the command with these inputs.</li>
   *   <li>Collects outputs in a `Map[String, Try[WdlValue]]`.</li>
   *   <li>If there are no `Failure`s among the outputs, messages the parent actor
   *   with `Completed`, otherwise messages `Failed`.</li>
   * </ol>
   */
  private def handleStart(): Unit = {
    // Record the original sender here and not in the yield over the `Future`, as the
    // value of sender() when that code executes may be different than what it is here.
    val originalSender = sender()
    val backendInputs = backend.adjustInputPaths(call, locallyQualifiedInputs)

    def handleFailedInstantiation(e: Throwable): Unit = {
      val message = s"Call '${call.fullyQualifiedName}' failed to launch command: " + e.getMessage
      log.error(s"$tag: $message")
      context.parent ! CallFailed(call, message)
    }

    def handleSuccessfulInstantiation(commandLine: String): Unit = {
      log.info(s"$tag: launching `$commandLine`")
      originalSender ! WorkflowActor.CallStarted(call)
      val tryOutputs = backend.executeCommand(commandLine, workflowDescriptor, call, s => locallyQualifiedInputs.get(s).get)
      val (successes, failures) = tryOutputs.partition {
        _._2.isSuccess
      }

      if (failures.isEmpty) {
        // Materialize the Successes.
        val outputs = successes.map { case (key, value) => key -> value.get }
        log.info(s"$tag: success")
        context.parent ! WorkflowActor.CallCompleted(call, outputs)
      } else {
        val errorMessages = TryUtil.stringifyFailures(failures.values).mkString("\n")
        log.error(s"$tag: failed")
        errorMessages foreach { m =>
          log.error(s"$tag: $m")
        }

        log.error(errorMessages.mkString("\n"))
        context.parent ! WorkflowActor.CallFailed(call, errorMessages.mkString("\n"))
      }
    }

    call.instantiateCommandLine(backendInputs) match {
      case Success(commandLine) => handleSuccessfulInstantiation(commandLine)
      case Failure(e) => handleFailedInstantiation(e)
    }
  }
}
