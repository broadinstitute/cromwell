package cromwell

import java.util.UUID

import cromwell.binding._
import cromwell.binding.types.WdlType
import cromwell.binding.values.WdlValue

/**
 * ==Cromwell Execution Engine==
 *
 * Given a WDL file and a backend to execute on, this package provides an API to launch a workflow
 * and query its progress.
 *
 * Internally, this package is built on top of [[cromwell.binding]].
 */
package object engine {
  case class WorkflowId(id: UUID) {
    override def toString = id.toString
    def shortString = id.toString.split("-")(0)
  }

  object WorkflowId {
    def fromString(id:String): WorkflowId = new WorkflowId(UUID.fromString(id))
    def randomId() = WorkflowId(UUID.randomUUID())
  }

  case class CallReference(workflowName: String, workflowId: WorkflowId, callName: String) {
    override def toString = s"UUID(${workflowId.shortString})/$callName"
  }

  case class AbortFunction(function: ()=>Unit)
  case class IndexedAbortFunction(callReference: CallReference, callAbortFunction: AbortFunction)
  case class AbortRegistrationFunction(register: AbortFunction=>Unit)

  sealed trait WorkflowState {
    def isTerminal: Boolean
  }

  private lazy val workflowStates = Seq(WorkflowSubmitted, WorkflowRunning, WorkflowFailed, WorkflowSucceeded, WorkflowAborting, WorkflowAborted)

  object WorkflowState {
    def fromString(str: String): WorkflowState = workflowStates.find(_.toString == str).getOrElse(
      throw new NoSuchElementException(s"No such WorkflowState: $str"))
  }

  case object WorkflowSubmitted extends WorkflowState {
    override val toString: String = "Submitted"
    override val isTerminal = false
  }
  
  case object WorkflowRunning extends WorkflowState {
    override val toString: String = "Running"
    override val isTerminal = false
  }

  case object WorkflowAborting extends WorkflowState {
    override val toString: String = "Aborting"
    override val isTerminal = false
  }

  case object WorkflowFailed extends WorkflowState {
    override val toString: String = "Failed"
    override val isTerminal = true
  }
  
  case object WorkflowSucceeded extends WorkflowState {
    override val toString: String = "Succeeded"
    override val isTerminal = true
  }

  case object WorkflowAborted extends WorkflowState {
    override val toString: String = "Aborted"
    override val isTerminal = true
  }

  case object WorkflowStopped extends WorkflowState {
    override val toString = "Stopped"
    override val isTerminal = true
  }

  object SymbolStoreEntry {
    private def splitFqn(fullyQualifiedName: FullyQualifiedName): (String, String) = {
      val lastIndex = fullyQualifiedName.lastIndexOf(".")
      (fullyQualifiedName.substring(0, lastIndex), fullyQualifiedName.substring(lastIndex + 1))
    }

    def apply(fullyQualifiedName: FullyQualifiedName, wdlValue: WdlValue, input: Boolean): SymbolStoreEntry = {
      val (scope, name) = splitFqn(fullyQualifiedName)
      val key = SymbolStoreKey(scope, name, index = None, input)
      SymbolStoreEntry(key, wdlValue.wdlType, Some(wdlValue))
    }

    def toWorkflowOutputs(t: Traversable[SymbolStoreEntry]): WorkflowOutputs = t.map { e =>
      s"${e.key.scope}.${e.key.name}" -> e.wdlValue.get
    }.toMap

    def toCallOutputs(traversable: Traversable[SymbolStoreEntry]): CallOutputs = traversable.map { entry =>
      entry.key.name -> entry.wdlValue.get
    }.toMap
  }

  case class SymbolStoreKey(scope: String, name: String, index: Option[Int], input: Boolean)

  case class SymbolStoreEntry(key: SymbolStoreKey, wdlType: WdlType, wdlValue: Option[WdlValue]) {
    def isInput: Boolean = key.input
    def isOutput: Boolean = !isInput
    def scope: String = key.scope
  }

  object ExecutionStatus extends Enumeration {
    type ExecutionStatus = Value
    val NotStarted, Starting, Running, Failed, Done, Aborted, Aborting = Value
  }
}
