package cromwell.engine

import cromwell.binding.{Call, WdlBinding}

object ExecutionStatus extends Enumeration {
  type ExecutionStatus = Value
  val NotStarted, Running, Failed, Done = Value
}


/**
 * Corresponds to the "execution table" of our discussions.
 */
class ExecutionStore(binding: WdlBinding) {

  def isWorkflowDone: Boolean = table.forall(_._2 == ExecutionStatus.Done)

  def updateStatus(call: Call, status: ExecutionStatus.Value): Unit = table += (call -> status)

  var table = binding.workflows.head.calls.map { call => call -> ExecutionStatus.NotStarted }.toMap

  /**
   * Return all calls which are currently in state `NotStarted` and whose prerequisites are all `Done`,
   * i.e. the calls which should now be eligible to run.
   */
  def runnableCalls: Iterable[Call] = {
    for {
      callEntry <- table if callEntry._2 == ExecutionStatus.NotStarted
      call = callEntry._1
      if call.prerequisiteCalls().forall {
        case c:Call => table.get(c).get == ExecutionStatus.Done
        case _ => false
      }
    } yield call
  }
  override def toString: String = {
    table.map {case(k, v) => s"${k.fullyQualifiedName}\t$v"}.mkString("\n")
  }
}
