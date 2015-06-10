package cromwell.engine

import akka.actor.{Actor, ActorRef, Props}
import akka.event.{Logging, LoggingReceive}
import akka.pattern.ask
import cromwell.binding.values.WdlValue
import cromwell.binding.{WorkflowDescriptor, Call}
import cromwell.engine.backend.Backend
import cromwell.util.TryUtil

import scala.concurrent.ExecutionContext.Implicits.global
import scala.concurrent.Future
import scala.language.postfixOps

object CallActor {
  sealed trait CallActorMessage
  case object Start extends CallActorMessage

  def props(call: Call, backend: Backend, workflowDescriptor: WorkflowDescriptor, storeActor: ActorRef): Props =
    Props(new CallActor(call, backend, workflowDescriptor, storeActor))
}


/** Actor to manage the execution of a single call. */
class CallActor(call: Call, backend: Backend, workflowDescriptor: WorkflowDescriptor, storeActor: ActorRef) extends Actor with CromwellActor {
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
   *   <li>Collects the locally qualified inputs from the store actor.</li>
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

    for {
      inputs <- (storeActor ? StoreActor.GetLocallyQualifiedInputs(call)).mapTo[Map[String, WdlValue]]
      backendInputs = backend.adjustInputPaths(call, inputs)
      commandLine <- Future.fromTry(call.instantiateCommandLine(backendInputs))
    } yield {
      log.info(s"$tag: launching `$commandLine`")
      originalSender ! WorkflowActor.CallStarted(call)
      val tryOutputs = backend.executeCommand(commandLine, workflowDescriptor, call, s => inputs.get(s).get)
      val (successes, failures) = tryOutputs.partition {
        _._2.isSuccess
      }

      if (failures.isEmpty) {
        // Materialize the Successes.
        val outputs = successes.map { case (key, value) => key -> value.get }
        log.info(s"$tag: success")
        context.parent ! WorkflowActor.CallCompleted(call, outputs)
      } else {
        val errorMessages = TryUtil.stringifyFailures(failures.values)
        log.error(s"$tag: failed")
        errorMessages foreach {m =>
          log.error(s"$tag: $m")
        }

        log.error(errorMessages.mkString("\n"))
        context.parent ! WorkflowActor.CallFailed(call, errorMessages.mkString("\n"))
      }
    }
  }
}
