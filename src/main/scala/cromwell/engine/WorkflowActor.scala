package cromwell.engine

import akka.actor.{Actor, ActorRef, Props}
import akka.event.{Logging, LoggingReceive}
import cromwell.binding.values.WdlValue
import cromwell.binding.{Call, FullyQualifiedName, WdlBinding}
import cromwell.engine.backend.Backend
import cromwell.parser.WdlParser.{Ast, Terminal}
import cromwell.util.TryUtil

import scala.concurrent.duration._
import scala.language.postfixOps


/**
 * This model assumes workflow actors are transient, created once per workflow instance.
 */
object WorkflowActor {

  sealed trait WorkflowActorMessage

  case class Start(backend: Backend) extends WorkflowActorMessage
  case object Started extends WorkflowActorMessage
  case class Done(symbolStore: SymbolStore) extends WorkflowActorMessage
  case class InvalidOperation(message: String) extends WorkflowActorMessage
  case class Failed(failure: String) extends WorkflowActorMessage
  case object CheckExecutionStatus extends WorkflowActorMessage

  def buildWorkflowActorProps(id: WorkflowId, binding: WdlBinding, actualInputs: Map[FullyQualifiedName, WdlValue]): Props = {
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

    Props(new WorkflowActor(id, binding, actualInputs))
  }
}

/** Represents the root of a single workflow instance, not a manager of multiple
  * workflows.
  */
case class WorkflowActor private(id: WorkflowId, binding: WdlBinding, actualInputs: Map[FullyQualifiedName, WdlValue]) extends Actor {

  import WorkflowActor._

  // That which started the workflow, a role which has not yet been defined more specifically.
  private var _primeMover: Option[ActorRef] = None

  override def receive = notStarted

  private var maybeSymbolStore: Option[SymbolStore] = None
  private var maybeExecutionStore: Option[ExecutionStore] = None
  private var maybeBackend: Option[Backend] = None
  private val log = Logging(context.system, this)


  private def notStarted: Receive = {
    LoggingReceive {

      case Start(backend) =>
        _primeMover = Option(sender())
        maybeExecutionStore = Option(new ExecutionStore(binding))
        maybeSymbolStore = Option(new SymbolStore(binding, actualInputs))
        maybeBackend = Option(backend)
        // Poll for runnable calls.
        context.system.scheduler.schedule(
          initialDelay = 0 milliseconds,
          interval = 5 seconds,
          receiver = self,
          message = CheckExecutionStatus)(context.system.dispatcher)

        context.become(started)
        primeMover ! Started


      case bad@_ =>
        sender ! InvalidOperation(s"Received message $bad before any Start message")
    }
  }

  private def createRequiredInputs(call: Call): Unit = {
    call.inputMappings.foreach { case(inputName, expression) =>
      val ast = expression.ast.asInstanceOf[Ast]
      val Seq(lhs, rhs) = Seq("lhs", "rhs").map { ast.getAttribute(_).asInstanceOf[Terminal].getSourceString }
      val outputFqn = Seq(call.parent.get.name, lhs, rhs).mkString(".")
      val inputFqn = Seq(call.parent.get.name, call.name, inputName).mkString(".")
      symbolStore.copyOutputToInput(outputFqn, inputFqn)
    }
  }

  private def started: Receive = {
    LoggingReceive {

      case CheckExecutionStatus =>
        if (executionStore.isDone) {
          primeMover ! Done(symbolStore)
        } else {
          val runnableCalls = executionStore.runnableCalls
          if (runnableCalls.nonEmpty) {
            val namesOfCallsBeingStarted = runnableCalls.map { _.name}.toSeq.sorted
            log.info("Starting calls: " + namesOfCallsBeingStarted.mkString(", "))
            runnableCalls.foreach { call =>
              createRequiredInputs(call)
              startCallActor(symbolStore, call)
            }
          }
        }

      case CallActor.Started(call) =>
        executionStore.updateStatus(call, ExecutionStatus.Running)

      case CallActor.Done(call, callOutputs) =>
        // TODO This should message store actors to synchronize updates.
        // TODO It doesn't seem right to write back any outputs if writing back some
        // TODO outputs fails.
        val addedEntries = callOutputs.map { callOutput =>
          symbolStore.addOutputValue(call.fullyQualifiedName, callOutput._1, Some(callOutput._2), callOutput._2.wdlType)
        }

        val failureMessages = TryUtil.stringifyFailures(addedEntries)
        if (failureMessages.isEmpty) {
          executionStore.updateStatus(call, ExecutionStatus.Done)
        } else {
          val errorMessages = failureMessages.mkString("\n")
          log.error(errorMessages)
          executionStore.updateStatus(call, ExecutionStatus.Failed)
        }
        self ! CheckExecutionStatus

      case CallActor.Failed(call, failure) =>
        executionStore.updateStatus(call, ExecutionStatus.Failed)
        primeMover ! Failed(failure)

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

  private def executionStore = maybeExecutionStore.get

  private def backend = maybeBackend.get
}
