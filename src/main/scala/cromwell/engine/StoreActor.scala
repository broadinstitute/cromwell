package cromwell.engine

import akka.actor.{Actor, Props}
import akka.event.{Logging, LoggingReceive}
import cromwell.binding._
import cromwell.binding.values.WdlValue
import cromwell.engine.ExecutionStatus.ExecutionStatus
import cromwell.engine.StoreActor._
import cromwell.util.TryUtil

import scala.language.postfixOps
import scala.util.Try

object StoreActor {
  def props(workflow: WorkflowDescriptor, hostInputs: HostInputs) = Props(new StoreActor(workflow, hostInputs))
  sealed trait StoreActorMessage
  case class CallCompleted(call: Call, callOutputs: Map[String, WdlValue]) extends StoreActorMessage
  case object StartRunnableCalls extends StoreActorMessage
  case class UpdateStatus(call: Call, status: ExecutionStatus.Value) extends StoreActorMessage
  case object GetOutputs extends StoreActorMessage
  case class GetLocallyQualifiedInputs(call: Call) extends StoreActorMessage
}

/**
 * Actor to hold symbol and execution status data for a single workflow.  This actor
 * guards mutable state over the symbol and execution stores. */
class StoreActor(workflow: WorkflowDescriptor, hostInputs: HostInputs) extends Actor with CromwellActor {
  private val symbolStore = new SymbolStore(workflow.namespace, hostInputs)
  private val executionStore = new ExecutionStore(workflow)
  private val log = Logging(context.system, this)
  val tag = s"StoreActor [UUID(${workflow.shortId})]"

  override def receive: Receive = LoggingReceive {
    case CallCompleted(call, callOutputs) =>
      updateOutputs(call, callOutputs)
      val msg = if (executionStore.isWorkflowDone) WorkflowActor.Complete else startRunnableCalls
      sender ! msg
    case StartRunnableCalls => sender ! startRunnableCalls
    case GetOutputs =>
      sender ! (symbolStore.getOutputs map symbolStoreEntryToMapEntry).toMap
    case GetLocallyQualifiedInputs(call) => sender ! symbolStore.locallyQualifiedInputs(call)
    case UpdateStatus(call, status) => executionStore.updateStatus(call, status)
  }

  private def symbolStoreEntryToMapEntry(e: SymbolStoreEntry): (String, WdlValue) = {
    e.key.scope + "." + e.key.name -> e.wdlValue.get
  }

  private def updateOutputs(call: Call, callOutputs: Map[String, WdlValue]): Unit = {
    def addOutputValueToSymbolStore(callOutput: (String, WdlValue)): Try[Unit] =
      symbolStore.addOutputValue(call.fullyQualifiedName, callOutput._1, Some(callOutput._2), callOutput._2.wdlType)

    callOutputs foreach {case (k, v) =>
      log.info(s"$tag: set ${call.fullyQualifiedName}.$k => $v")
    }

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

  private def startRunnableCalls: WorkflowActor.RunnableCalls = WorkflowActor.RunnableCalls(executionStore.startRunnableCalls)
}
