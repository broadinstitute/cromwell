package cromwell.engine

import akka.actor.{Actor, Props}
import akka.event.{Logging, LoggingReceive}
import cromwell.binding.values.WdlValue
import cromwell.binding.{Call, CallInputs, WorkflowDescriptor}
import cromwell.engine.backend.Backend
import cromwell.engine.workflow.WorkflowActor
import cromwell.util.TryUtil

import scala.language.postfixOps


object CallActor {
  sealed trait CallActorMessage
  case class Start(inputs: CallInputs) extends CallActorMessage

  def props(call: Call, locallyQualifiedInputs: Map[String, WdlValue], backend: Backend, workflowDescriptor: WorkflowDescriptor, name: String): Props =
    Props(new CallActor(call, locallyQualifiedInputs, backend, workflowDescriptor))
}


/** Actor to manage the execution of a single call. */
class CallActor(call: Call, locallyQualifiedInputs: Map[String, WdlValue], backend: Backend, workflowDescriptor: WorkflowDescriptor) extends Actor with CromwellActor {

  private val log = Logging(context.system, this)

  override def receive = LoggingReceive {
    case CallActor.Start => handleStart()
    case badMessage =>
      val diagnostic = s"Received unexpected message $badMessage."
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
    log.info("In handleStart for " + call.fullyQualifiedName)
    val originalSender = sender()
    val backendInputs = backend.adjustInputPaths(call, locallyQualifiedInputs)
    log.info("Got backend inputs for " + call.fullyQualifiedName + ": " + backendInputs)

    val tryCommand = for {
      commandLine <- call.instantiateCommandLine(backendInputs)
    } yield {
      log.info("Command line for " + call.fullyQualifiedName + ": " + commandLine)
      originalSender ! WorkflowActor.CallStarted(call)
      val tryOutputs = backend.executeCommand(commandLine, workflowDescriptor, call, s => locallyQualifiedInputs.get(s).get)
      log.info(s"Executed command '$commandLine'")
      val (successes, failures) = tryOutputs.partition {
        _._2.isSuccess
      }

      if (failures.isEmpty) {
        // Materialize the Successes.
        val outputs = successes.map { case (key, value) => key -> value.get }
        context.parent ! WorkflowActor.CallCompleted(call, outputs)
      } else {
        val errorMessages = TryUtil.stringifyFailures(failures.values).mkString("\n")

        log.error(errorMessages)
        context.parent ! WorkflowActor.CallFailed(call, errorMessages)
      }
    }
    tryCommand.recover {
      case e: Throwable => log.error(e, s"Call '${call.fullyQualifiedName}' failed to execute command: " + e.getMessage)
    }
  }
}
