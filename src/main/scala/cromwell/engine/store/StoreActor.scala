package cromwell.engine.store

import scala.language.postfixOps
import scala.util.Try
import akka.actor.{Actor, Props}
import akka.event.{Logging, LoggingReceive}
import cromwell.binding._
import cromwell.binding.values.WdlValue
import cromwell.engine._
import cromwell.engine.db.CallInfo
import cromwell.engine.store.ExecutionStore.ExecutionStatus
import cromwell.engine.store.StoreActor.InitialStore
import cromwell.engine.workflow.WorkflowActor
import cromwell.util.TryUtil
import StoreActor._

object StoreActor {
  case class InitialStore(symbolStore: Set[SymbolStoreEntry], executionStore: Set[CallInfo])

  /** Creates a StoreActor from previous data */
  def props(namespace: WdlNamespace, initialStore: InitialStore): Props = {
    Props(new StoreActor(namespace, SymbolStore(initialStore.symbolStore), ExecutionStore(namespace, initialStore.executionStore)))
  }

  /** Constructs a StoreActor from scratch */
  def props(namespace: WdlNamespace, hostInputs: HostInputs): Props = {
    Props(new StoreActor(namespace, SymbolStore(namespace, hostInputs), ExecutionStore(namespace)))
  }

  /**
   * If initialStore is defined uses that otherwise hostInputs. Never uses information from both arguments.
   *
   * FIXME: Perhaps the two args should be an Either?
   */
  def props(namespace: WdlNamespace, hostInputs: HostInputs, initialStore: Option[InitialStore]): Props =  {
    initialStore map {props(namespace, _)} getOrElse props(namespace, hostInputs)
  }

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
 class StoreActor(namespace: WdlNamespace,
                  private val symbolStore: SymbolStore,
                  private val executionStore: ExecutionStore) extends Actor with CromwellActor {
  private val log = Logging(context.system, this)

  override def receive: Receive = LoggingReceive {
    case CallCompleted(call, callOutputs) =>
      updateOutputs(call, callOutputs)
      val msg = if (executionStore.isWorkflowDone) WorkflowActor.Complete else startRunnableCalls
      sender ! msg
    case StartRunnableCalls => sender ! startRunnableCalls
    case GetOutputs =>
      sender ! (symbolStore.outputs map symbolStoreEntryToMapEntry).toMap
    case GetLocallyQualifiedInputs(call) => sender ! symbolStore.locallyQualifiedInputs(call)
    case UpdateStatus(call, status) => executionStore.updateStatus(call, status)
  }

  private def symbolStoreEntryToMapEntry(e: SymbolStoreEntry): (String, WdlValue) = {
    e.key.scope + "." + e.key.name -> e.wdlValue.get
  }
  
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

  private def startRunnableCalls = WorkflowActor.RunnableCalls(executionStore.startRunnableCalls)
}
