package cromwell.backend

import cromwell.backend.BackendExecutionStatus.BackendExecutionStatus

// TODO: PBE: Not moving the "engine" version yet as it has many other "engine" classes mixed in. Eventually "core"?

object BackendExecutionStore {
  def empty = BackendExecutionStore(Map.empty)
  type ExecutionStoreEntry = (JobKey, BackendExecutionStatus)
}

case class BackendExecutionStore(store: Map[JobKey, BackendExecutionStatus])
