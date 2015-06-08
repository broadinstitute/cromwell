package cromwell.engine

import akka.actor.{Actor, Props}
import akka.event.{Logging, LoggingReceive}
import akka.util.Timeout
import cromwell.binding.values.WdlValue
import cromwell.binding.{WdlBinding, _}
import cromwell.engine.StoreActor._
import cromwell.parser.WdlParser.{Ast, Terminal}
import cromwell.util.TryUtil

import scala.concurrent.duration._
import scala.language.postfixOps
import scala.util.Try


object StoreActor {
  def props(binding: WdlBinding, inputs: WorkflowCoercedInputs) =
    Props(new StoreActor(binding, inputs))

  sealed trait StoreActorMessage
  case class CallCompleted(call: Call, callOutputs: Map[String, WdlValue]) extends StoreActorMessage
  case object FindRunnableCalls extends StoreActorMessage
  case class UpdateStatus(call: Call, status: ExecutionStatus.Value) extends StoreActorMessage
  case object GetOutputs extends StoreActorMessage
  case class GetLocallyQualifiedInputs(call: Call) extends StoreActorMessage
  implicit val ActorTimeout = Timeout(5 seconds)
}

/**
 * Actor to hold symbol and execution status data for a single workflow.  This actor
 * guards mutable state over the symbol and execution stores, and must therefore not
 * pass back `Future`s over updates or reads of those stores. */
class StoreActor(binding: WdlBinding, inputs: WorkflowCoercedInputs) extends Actor {
  private val symbolStore = new SymbolStore(binding, inputs)
  private val executionStore = new ExecutionStore(binding)
  private val log = Logging(context.system, this)

  override def receive: Receive = LoggingReceive {

    case CallCompleted(call, callOutputs) =>
      sender ! handleCallCompleted(call, callOutputs)

    case FindRunnableCalls =>
      sender ! executionStore.runnableCalls

    case GetOutputs =>
      sender ! (symbolStore.getOutputs map symbolStoreEntryToMapEntry).toMap

    case GetLocallyQualifiedInputs(call) =>
      sender ! symbolStore.locallyQualifiedInputs(call)

    case UpdateStatus(call, status) =>
      executionStore.updateStatus(call, status)
  }

  private def symbolStoreEntryToMapEntry(e: SymbolStoreEntry): (String, WdlValue) =
    e.key.scope + "." + e.key.name -> e.wdlValue.get

  private def updateOutputs(call: Call, callOutputs: Map[String, WdlValue]): Unit = {

    def addOutputValueToSymbolStore(callOutput: (String, WdlValue)): Try[Unit] =
      symbolStore.addOutputValue(call.fullyQualifiedName, callOutput._1, Some(callOutput._2), callOutput._2.wdlType)

    val addedEntries = callOutputs map addOutputValueToSymbolStore
    val failureMessages = TryUtil.stringifyFailures(addedEntries)

    if (failureMessages.isEmpty) {
      executionStore.updateStatus(call, ExecutionStatus.Done)
    } else {
      val errorMessages = failureMessages.mkString("\n")
      log.error(errorMessages)
      executionStore.updateStatus(call, ExecutionStatus.Failed)
      throw new RuntimeException(errorMessages)
    }
  }

  /**
   * Updates outputs for the completed call and returns a `Future[Boolean]` which is true
   * if the workflow is now "done".  Current "done" for a workflow means all calls are done.
   */
  private def handleCallCompleted(call: Call, callOutputs: Map[String, WdlValue]): Boolean = {
    updateOutputs(call, callOutputs)
    executionStore.isWorkflowDone
  }

}
