package cromwell.engine

import cromwell.engine.ExecutionIndex.ExecutionIndex
import cromwell.engine.ExecutionStatus.ExecutionStatus
import cromwell.engine.backend.{ReturnCode, BackendCall}

package object db {
  case class CallStatus(executionStatus: ExecutionStatus, returnCode: Option[ReturnCode], hash: Option[ExecutionHash], resultsClonedFrom: Option[BackendCall]) {
    def isTerminal: Boolean = executionStatus.isTerminal
    def isStarting: Boolean = executionStatus == ExecutionStatus.Starting
  }

  // Uniquely identify an entry in the execution table
  case class ExecutionDatabaseKey(fqn: FullyQualifiedName, index: ExecutionIndex, attempt: Int = 1) {
    def isCollector(keys: Traversable[ExecutionDatabaseKey]): Boolean = {
      index.isEmpty &&
        (keys exists { e =>
          (e.fqn == fqn) && e.index.isDefined
        })
    }
  }

  case class JesId(id: String)
  case class JesStatus(status: String)
}
