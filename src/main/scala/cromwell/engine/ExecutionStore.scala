package cromwell.engine

import cromwell.binding.{Call, WdlNamespace}
import cromwell.engine.db.CallInfo

object ExecutionStatus extends Enumeration {
  type ExecutionStatus = Value
  val NotStarted, Starting, Running, Failed, Done = Value

  def apply(namespace: WdlNamespace): ExecutionStore = new ExecutionStore(namespace, Set.empty[CallInfo])
}

/**
 * Corresponds to the "execution table" of our discussions.
 */
case class ExecutionStore(namespace: WdlNamespace, initialStore: Set[CallInfo]) {
  // FIXME: this makes 'calls' passes over the initial store, perhaps partition and do things in one shot?
  private var table = namespace.workflows.head.calls.map {determineCallStatus}.toMap

  /**
   * Start all calls which are currently in state `NotStarted` and whose prerequisites are all `Done`,
   * i.e. the calls which should now be eligible to run.
   */
  def startRunnableCalls: Iterable[Call] = {
    val runnableCalls = for {
      callEntry <- table if callEntry._2 == ExecutionStatus.NotStarted
      call = callEntry._1 if callIsRunnable(call)
    } yield call

    runnableCalls foreach { call =>
      table += call -> ExecutionStatus.Starting
    }
    runnableCalls
  }

  def isWorkflowDone: Boolean = table.forall(_._2 == ExecutionStatus.Done)
  def updateStatus(call: Call, status: ExecutionStatus.Value): Unit = table += (call -> status)

  private def callIsRunnable(call: Call): Boolean = {
    call.prerequisiteCalls().forall(table.get(_).get == ExecutionStatus.Done)
  }

  private def determineCallStatus(call: Call) = {
    val status = initialStore.find(_.callFqn == call.taskFqn) map {_.status} getOrElse ExecutionStatus.NotStarted
    call -> status
  }

  override def toString: String = {
    table.map {case(k, v) => s"${k.fullyQualifiedName}\t$v"}.mkString("\n")
  }
}
