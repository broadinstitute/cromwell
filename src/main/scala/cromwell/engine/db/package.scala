package cromwell.engine

import cromwell.binding.FullyQualifiedName
import cromwell.engine.ExecutionIndex.ExecutionIndex
import cromwell.engine.ExecutionStatus.ExecutionStatus

package object db {
  // CallStatus doesn't exist as an enum already, or reuse ExecutionStatus?
  type CallStatus = ExecutionStatus
  //Uniquely identify an entry in the execution table
  case class ExecutionDatabaseKey(fqn: FullyQualifiedName, index: ExecutionIndex)
  // jesId and jesStatus should have stronger types
  type JesId = Int
  type JesStatus = String
}
