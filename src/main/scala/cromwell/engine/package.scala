package cromwell

import java.util.UUID

import cromwell.binding.{WdlNamespace, WorkflowCoercedInputs}

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
}
