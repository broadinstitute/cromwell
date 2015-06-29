package cromwell.engine

import akka.event.Logging
import cromwell.binding.{WorkflowDescriptor, Call, WdlNamespace}
import org.slf4j.LoggerFactory

object ExecutionStatus extends Enumeration {
  type ExecutionStatus = Value
  val NotStarted, Starting, Running, Failed, Done = Value
}


/**
 * Corresponds to the "execution table" of our discussions.
 */
class ExecutionStore(workflow: WorkflowDescriptor) {
  val log = LoggerFactory.getLogger("ExecutionStore")
  val tag = s"ExecutionStore [UUID(${workflow.shortId})]"
  var table = workflow.namespace.workflows.head.calls.map { call => call -> ExecutionStatus.NotStarted }.toMap
  def isWorkflowDone: Boolean = table.forall(_._2 == ExecutionStatus.Done)

  def updateStatus(call: Call, status: ExecutionStatus.Value): Unit = {
    log.info(s"$tag: ${call.name} update to $status")
    table += (call -> status)
  }

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

  private def callIsRunnable(call: Call): Boolean = {
    call.prerequisiteCalls().forall(table.get(_).get == ExecutionStatus.Done)
  }

  override def toString: String = {
    table.map {case(k, v) => s"${k.fullyQualifiedName}\t$v"}.mkString("\n")
  }
}
