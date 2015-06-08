package cromwell.engine

import akka.actor.{Actor, Props}
import akka.event.{Logging, LoggingReceive}
import akka.pattern.{ask, pipe}
import akka.util.Timeout
import cromwell.binding.values.WdlValue
import cromwell.binding.{Call, WdlBinding, WorkflowCoercedInputs, WorkflowOutputs}
import cromwell.engine.StoreActor._
import cromwell.engine.backend.Backend

import scala.concurrent.ExecutionContext.Implicits.global
import scala.concurrent.duration._
import scala.language.postfixOps
import scala.util.{Failure, Success}


object WorkflowActor {

  sealed trait WorkflowActorMessage

  case object Start extends WorkflowActorMessage
  case object Started extends WorkflowActorMessage
  case class Done(outputs: Map[String, WdlValue]) extends WorkflowActorMessage
  case class Failed(failure: String) extends WorkflowActorMessage
  case object GetOutputs extends WorkflowActorMessage

  def props(id: WorkflowId, binding: WdlBinding, coercedInputs: WorkflowCoercedInputs, backend: Backend) =
    Props(new WorkflowActor(id, binding, coercedInputs, backend))

  implicit val ActorTimeout = Timeout(5 seconds)
}

/**
 * Actor to manage the execution of a single workflow.
 */
case class WorkflowActor private(id: WorkflowId,
                                 binding: WdlBinding,
                                 actualInputs: WorkflowCoercedInputs,
                                 backend: Backend) extends Actor {

  private val log = Logging(context.system, this)

  private val storeActor = context.actorOf(StoreActor.props(binding, actualInputs))

  def receive: Receive = unstarted

  /** Unstarted `Receive` state that only responds to `Start` messages. */
  def unstarted: Receive = LoggingReceive {
    case WorkflowActor.Start => handleInitialStart()
    case bad => fail(bad, "unstarted")
  }

  /** Started `Receive` state. */
  def started: Receive = LoggingReceive {
    case CallActor.Started(call) =>
      storeActor ! UpdateStatus(call, ExecutionStatus.Running)

    case CallActor.Completed(completedCall, callOutputs) =>
      handleCallCompleted(completedCall, callOutputs)

    case CallActor.Failed(call, failure) =>
      storeActor ! UpdateStatus(call, ExecutionStatus.Failed)
      context.parent ! WorkflowActor.Failed(failure)

    case WorkflowActor.GetOutputs =>
      (storeActor ? StoreActor.GetOutputs) pipeTo sender

    case bad => fail(bad, "started")
  }

  /** Message the parent actor with `Failed` with an appropriate message. */
  private def fail(badMessage: Any, stateName: String): Unit = {
    val diagnostic = s"Unexpected message $badMessage received when in $stateName state."
    log.error(diagnostic)
    context.parent ! WorkflowActor.Failed(diagnostic)
  }

  /**
   * Performs the following steps:
   * <ol>
   *   <li>Finds all runnable calls.</li>
   *   <li>Copies outputs from prerequisite calls to the inputs of the runnable calls.</li>
   *   <li>Starts a call actor for each runnable call.</li>
   * </ol>
   */
  private def startRunnableCalls(): Unit = {
    for {
      runnableCalls <- (storeActor ? FindRunnableCalls).mapTo[Iterable[Call]]
    } yield {
      if (runnableCalls.nonEmpty) {
        log.info("Starting calls: " + runnableCalls.map {_.name}.toSeq.sorted.mkString(", "))
      }
      runnableCalls foreach startCallActor
    }
  }

  /**
   * Performs the following steps:
   * <ol>
   *   <li>Puts this actor into its "started" state.</li>
   *   <li>Starts all runnable calls.</li>
   *   <li>Messages the sender that the workflow has been started.</li>
   * <ol>
   */
  private def handleInitialStart(): Unit = {
    context.become(started)
    startRunnableCalls()
    sender ! WorkflowActor.Started
  }

  /**
   * Coordinate handling of the completion of this call with the store actor.  The invocation of this
   * method should result in the outputs of the call being written to the symbol store and the status
   * of the call being marked complete.  If all calls in the workflow are complete the parent of this
   * actor will be messaged a `WorkflowActor.Done` with the workflow outputs.
   */
  private def handleCallCompleted(completedCall: Call, callOutputs: WorkflowOutputs): Unit = {
    // Messages the store actor with a `CallCompleted` via an `ask`; i.e. this is not *asking* the
    // store actor if the call is completed, it's informing the store actor of that fact via `ask`.
    // The store actor responds with a `Future[Boolean]` which if `true` means the workflow has completed.
    (storeActor ? CallCompleted(completedCall, callOutputs)) onComplete {
      case Success(done) =>
        done match {
          case true =>
            val futureDoneMessage = for {
              outputs <- (storeActor ? GetOutputs).mapTo[Map[String, WdlValue]]
            } yield WorkflowActor.Done(outputs)
            futureDoneMessage pipeTo context.parent

          case _ => startRunnableCalls()
        }
      case Failure(t) =>
        log.error(t, t.getMessage)
        context.parent ! WorkflowFailed
    }
  }

  /** Create a per-call `CallActor` for the specified `Call` and send it a `Start` message to
    * begin execution. */
  private def startCallActor(call: Call): Unit = {
    val callActorProps = CallActor.props(call, backend, storeActor, "CallActor-" + call.name)
    context.actorOf(callActorProps) ! CallActor.Start
  }
}
