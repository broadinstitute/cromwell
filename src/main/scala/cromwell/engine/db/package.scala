package cromwell.engine

import cromwell.engine.store.ExecutionStore
import ExecutionStore.ExecutionStatus.ExecutionStatus

package object db {
  // CallStatus doesn't exist as an enum already, or reuse ExecutionStatus?
  type CallStatus = ExecutionStatus
  // jesId and jesStatus should have stronger types
  type JesId = String
  type JesStatus = String
}
