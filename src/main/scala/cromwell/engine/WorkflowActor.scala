package cromwell.engine

import akka.actor.SupervisorStrategy.Stop
import akka.actor.{Actor, ActorRef, Props}
import akka.event.{Logging, LoggingReceive}
import cromwell.binding.values.WdlValue
import cromwell.binding.{Call, FullyQualifiedName, WdlBinding}
import cromwell.engine.backend.Backend
import cromwell.util.TryUtil

/**
 * This model assumes workflow actors are transient, created once per workflow instance.
 */
object WorkflowActor {

  sealed trait WorkflowActorMessage

  case class Start(backend: Backend) extends WorkflowActorMessage
  case object Started extends WorkflowActorMessage

  // Note there is no Stop message defined here, Stop is from Akka.
  case object Stopped extends WorkflowActorMessage
  case class Done(symbolStore: SymbolStore) extends WorkflowActorMessage
  case class InvalidOperation(message: String) extends WorkflowActorMessage
  case class Failed(failure: String) extends WorkflowActorMessage

  def buildWorkflowActorProps(binding: WdlBinding, actualInputs: Map[FullyQualifiedName, WdlValue]): Props = {
    // Check inputs for presence and type compatibility
    val diagnostics = binding.workflow.inputs.collect {
      case requiredInput if !actualInputs.contains(requiredInput._1) =>
        requiredInput._1 -> "Required workflow input not specified"

      case requiredInput if actualInputs.get(requiredInput._1).get.wdlType != requiredInput._2 =>
        val expected = actualInputs.get(requiredInput._1).get.wdlType
        // FIXME formatting
        requiredInput._1 -> s"Incompatible workflow input types, expected $expected, got ${requiredInput._2}"
    }

    if (diagnostics.nonEmpty) throw new UnsatisfiedInputsException(diagnostics)

    Props(new WorkflowActor(binding, actualInputs))
  }
}

/** Represents the root of a single workflow instance, not a manager of multiple
  * workflows.
  */
class WorkflowActor private(binding: WdlBinding, actualInputs: Map[FullyQualifiedName, WdlValue]) extends Actor {

  import WorkflowActor._

  // That which started the workflow, a role which has not yet been defined more specifically.
  private var _primeMover: Option[ActorRef] = None

  override def receive = notStarted

  private var maybeSymbolStore: Option[SymbolStore] = None
  private var maybeExecutionStatusStore: Option[ExecutionStatusStore] = None
  private var maybeBackend: Option[Backend] = None
  private val log = Logging(context.system, this)


  private def notStarted: Receive = {
    LoggingReceive {

      case Start(backend) =>
        _primeMover = Option(sender())
        maybeExecutionStatusStore = Option(new ExecutionStatusStore(binding))
        maybeSymbolStore = Option(new SymbolStore(binding, actualInputs))
        maybeBackend = Option(backend)
        executionStatusStore.runnableCalls.foreach { call =>
          startCallActor(symbolStore, call)
        }
        primeMover ! Started
        context.become(started)

      case bad @ _ =>
        sender ! InvalidOperation(s"Received message $bad before any Start message")
    }
  }

  private def started: Receive = {
    LoggingReceive {
      case CallActor.Started =>
      // TODO update execution status store

      case CallActor.Done(call, callOutputs) =>
        // TODO Update execution status and symbol stores.
        // TODO There is only one runnable task at the moment, so immediately
        // TODO report as done to the creator, passing back the symbol store.
        val addedEntries = callOutputs.map { callOutput =>
          // TODO very broken, assumes there are always output values here, but other code
          // TODO has them properly Option'ed.
          // TODO This should be messaging a symbol store actor to synchronize updates
          symbolStore.addOutputValue(call.fullyQualifiedName, callOutput._1, Some(callOutput._2), callOutput._2.wdlType)
        }

        val failureMessages = TryUtil.stringifyFailures(addedEntries)
        if (failureMessages.isEmpty) {
          primeMover ! Done(symbolStore)
        } else {
          val errorMessages = failureMessages.mkString("\n")
          log.error(errorMessages)
          primeMover ! Failed(errorMessages)
        }
      // TODO This should reevaluate runnable calls, if there are no more
      // TODO initiate shutdown, including waiting for child Call actors to
      // TODO report as Stopped (state "SHUTTING_DOWN"?).  Shutting this down
      // TODO before waiting for child actors to report as stopped
      // TODO would produce a slew of dead letter warnings as a shutdown
      // TODO of self would terminate the entire actor hierarchy.
      // TODO If there are more runnable tasks they should be started and
      // TODO the stores updated appropriately.
      // TODO The code needed here is not very different from what Start
      // TODO should be doing.

      case CallActor.Failed(failure) =>
        primeMover ! Failed(failure)

      case CallActor.Stopped =>
      // TODO update execution status store

      case Stop =>
        // TODO is this the right way to stop this actor?
        // TODO http://stackoverflow.com/questions/13847963/akka-kill-vs-stop-vs-poison-pill
        context.stop(self)
        primeMover ! Stopped

      case unknown@_ =>
        primeMover ! InvalidOperation(s"Unknown/invalid message '$unknown'")
    }
  }

  private def primeMover = _primeMover.get

  private def startCallActor(symbolStore: SymbolStore, call: Call): Unit = {
    val callActor = context.actorOf(CallActor.props, "CallActor-" + call.name)
    callActor ! CallActor.Start(call, backend, symbolStore)
  }

  private def symbolStore = maybeSymbolStore.get

  private def executionStatusStore = maybeExecutionStatusStore.get

  private def backend = maybeBackend.get
}
