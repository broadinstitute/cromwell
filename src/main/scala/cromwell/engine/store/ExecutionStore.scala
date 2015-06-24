package cromwell.engine.store

import cromwell.binding.{Call, WdlNamespace}
import cromwell.engine.db.CallInfo
import cromwell.engine.store.ExecutionStore.ExecutionStatus

object ExecutionStore {
  object ExecutionStatus extends Enumeration {
    type ExecutionStatus = Value
    val NotStarted, Starting, Running, Failed, Done = Value
  }

  def apply(namespace: WdlNamespace): ExecutionStore = new ExecutionStore(namespace, Set.empty[CallInfo])
}

/** Contains a mapping of calls to their current status */
case class ExecutionStore(namespace: WdlNamespace, initialStore: Set[CallInfo]) {
  private var store: Map[Call, ExecutionStatus.Value] = createStore()

  /**
   * Start all calls which are currently in state `NotStarted` and whose prerequisites are all `Done`,
   * i.e. the calls which should now be eligible to run.
   */
  def startRunnableCalls: Iterable[Call] = {
    val runnableCalls = for {
      callEntry <- store if callEntry._2 == ExecutionStatus.NotStarted
      call = callEntry._1 if callIsRunnable(call)
    } yield call

    runnableCalls foreach {c => updateStatus(c, ExecutionStatus.Starting)}

    runnableCalls
  }

  def isWorkflowDone: Boolean = store.forall(_._2 == ExecutionStatus.Done)

  def updateStatus(call: Call, status: ExecutionStatus.Value): Unit = store += (call -> status)

  private def callIsRunnable(call: Call) = call.prerequisiteCalls().forall(store.get(_).get == ExecutionStatus.Done)

  private def createStore(): Map[Call, ExecutionStatus.Value] = {
    val nameToStatus = initialStore.map(c => c.callFqn -> c.status).toMap
    namespace.workflows.head.calls.map {c => c -> nameToStatus.getOrElse(c.taskFqn, ExecutionStatus.NotStarted)}.toMap
  }

  override def toString: String = store.map {case(k, v) => s"${k.fullyQualifiedName}\t$v"}.mkString("\n")
}
