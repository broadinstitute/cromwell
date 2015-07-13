package cromwell.engine

import cromwell.engine.ExecutionStatus.ExecutionStatus

package object db {
  // CallStatus doesn't exist as an enum already, or reuse ExecutionStatus?
  type CallStatus = ExecutionStatus
  // jesId and jesStatus should have stronger types
  type JesId = Int
  type JesStatus = String
}
