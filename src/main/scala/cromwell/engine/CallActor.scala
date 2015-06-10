package cromwell.engine

import akka.actor.{Actor, ActorRef, Props}
import akka.event.{Logging, LoggingReceive}
import akka.pattern.ask
import akka.util.Timeout
import cromwell.binding.values.WdlValue
import cromwell.binding.{Call, WorkflowOutputs}
import cromwell.engine.backend.Backend
import cromwell.util.TryUtil

import scala.concurrent.ExecutionContext.Implicits.global
import scala.concurrent.Future
import scala.concurrent.duration._
import scala.language.postfixOps


object CallActor {
  sealed trait CallActorMessage
  case object Start extends CallActorMessage
  case class Started(call: Call) extends CallActorMessage
  case class Completed(call: Call, outputs: WorkflowOutputs) extends CallActorMessage
  case class Failed(call: Call, failure: String) extends CallActorMessage

  def props(call: Call, backend: Backend, storeActor: ActorRef, name: String): Props =
    Props(new CallActor(call, backend, storeActor))

  implicit val ActorTimeout = Timeout(5 seconds)
}


/** Actor to manage the execution of a single call. */
class CallActor(call: Call, backend: Backend, storeActor: ActorRef) extends Actor {
  import CallActor.ActorTimeout

  private val log = Logging(context.system, this)

  override def receive = LoggingReceive {
    case CallActor.Start => handleStart()
    case badMessage =>
      val diagnostic = s"Received unexpected message $badMessage."
      log.error(diagnostic)
      context.parent ! CallActor.Failed(call, diagnostic)
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
      commandLine <- Future.fromTry(call.task.command.instantiate(inputs))
    } yield {
      originalSender ! CallActor.Started(call)
      val tryOutputs = backend.executeCommand(commandLine, call, call.task.outputs, s => inputs.get(s).get)
      val (successes, failures) = tryOutputs.partition {
        _._2.isSuccess
      }

      if (failures.isEmpty) {
        // Materialize the Successes.
        val outputs = successes.map { case (key, value) => key -> value.get }
        context.parent ! CallActor.Completed(call, outputs)
      } else {
        val errorMessages = TryUtil.stringifyFailures(failures.values).mkString("\n")

        log.error(errorMessages)
        context.parent ! CallActor.Failed(call, errorMessages)
      }
    }
  }
}
