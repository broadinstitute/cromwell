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
  case class ExecutionDatabaseKey(fqn: FullyQualifiedName, index: ExecutionIndex)
  // jesId and jesStatus should have stronger types
  type JesId = Int
  type JesStatus = String
}
