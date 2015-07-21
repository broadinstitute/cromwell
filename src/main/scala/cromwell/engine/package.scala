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
  type WorkflowId = UUID

  sealed trait WorkflowState {
    def isTerminal: Boolean
  }

  private lazy val workflowStates = Seq(WorkflowSubmitted, WorkflowRunning, WorkflowFailed, WorkflowSucceeded)

  object WorkflowState {
    def fromString(str: String): WorkflowState = workflowStates.find(_.toString == str).getOrElse(
      throw new NoSuchElementException(s"No such WorkflowState: $str"))
  }

  case object WorkflowSubmitted extends WorkflowState {
    override def toString: String = "Submitted"
    override val isTerminal = false
  }
  case object WorkflowRunning extends WorkflowState {
    override def toString: String = "Running"
    override val isTerminal = false
  }
  case object WorkflowFailed extends WorkflowState {
    override def toString: String = "Failed"
    override val isTerminal = true
  }
  case object WorkflowSucceeded extends WorkflowState {
    override def toString: String = "Succeeded"
    override val isTerminal = true
  }

  object ExecutionStatus extends Enumeration {
    type ExecutionStatus = Value
    val NotStarted, Starting, Running, Failed, Done = Value
  }

  object SymbolStoreEntry {
    private def splitFqn(fullyQualifiedName: FullyQualifiedName): (String, String) = {
      val lastIndex = fullyQualifiedName.lastIndexOf(".")
      (fullyQualifiedName.substring(0, lastIndex), fullyQualifiedName.substring(lastIndex + 1))
    }

    def apply(fullyQualifiedName: FullyQualifiedName, wdlValue: WdlValue, input: Boolean): SymbolStoreEntry = {
      val (scope, name) = splitFqn(fullyQualifiedName)
      val key = SymbolStoreKey(scope, name, iteration = None, input)
      SymbolStoreEntry(key, wdlValue.wdlType, Some(wdlValue))
    }

    def toWorkflowOutputs(t: Traversable[SymbolStoreEntry]): WorkflowOutputs = t.map { e =>
      s"${e.key.scope}.${e.key.name}" -> e.wdlValue.get
    }.toMap
  }

  case class SymbolStoreKey(scope: String, name: String, iteration: Option[Int], input: Boolean)

  case class SymbolStoreEntry(key: SymbolStoreKey, wdlType: WdlType, wdlValue: Option[WdlValue]) {
    def isInput: Boolean = key.input
    def isOutput: Boolean = !isInput
    def scope: String = key.scope
  }
}
