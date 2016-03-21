package cromwell.engine

import cromwell.engine.ExecutionStatus.ExecutionStatus
import cromwell.engine.backend.BackendCallJobDescriptor

package object db {
  case class CallStatus(executionStatus: ExecutionStatus, returnCode: Option[Int], hash: Option[ExecutionHash], resultsClonedFrom: Option[BackendCallJobDescriptor]) {
    def isTerminal: Boolean = executionStatus.isTerminal
    def isStarting: Boolean = executionStatus == ExecutionStatus.Starting
  }

  case class JesId(id: String)
  case class JesStatus(status: String)
}
