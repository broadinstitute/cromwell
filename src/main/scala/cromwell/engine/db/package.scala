package cromwell.engine

import cromwell.binding.FullyQualifiedName
import cromwell.engine.ExecutionIndex.ExecutionIndex
import cromwell.engine.ExecutionStatus.ExecutionStatus

package object db {
  case class CallStatus(executionStatus: ExecutionStatus, rc: Option[Int])

  object CallStatus {
    def apply(statusName: String, rc: Option[Int] = None): CallStatus = CallStatus(ExecutionStatus.withName(statusName), rc)
  }

  // Uniquely identify an entry in the execution table
  case class ExecutionDatabaseKey(fqn: FullyQualifiedName, index: ExecutionIndex) {
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
