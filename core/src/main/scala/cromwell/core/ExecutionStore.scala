package cromwell.core

import cromwell.core.ExecutionStatus._


object ExecutionStore {
  def empty = ExecutionStore(Map.empty)
  type ExecutionStoreEntry = (JobKey, ExecutionStatus)
}

case class ExecutionStore(store: Map[JobKey, ExecutionStatus]) {
  def add(values: Map[JobKey, ExecutionStatus]) = this.copy(store = store ++ values)
}
