package cromwell.engine

import cromwell.backend.model.ExecutionHash
import cromwell.engine.ExecutionIndex.ExecutionIndex
import cromwell.engine.ExecutionStatus.ExecutionStatus

package object db {
  case class CallStatus(executionStatus: ExecutionStatus, returnCode: Option[Int], hash: Option[ExecutionHash]) {
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
