package cromwell.backend.model

import cromwell.backend.model.Status.Status

/**
  * Represents status of an execution.
  */
object Status extends Enumeration {
  type Status = Value
  val Created, Running, Succeeded, Failed, Canceled = Value
}

/**
  * Defines a task status with its result data.
  */
case class TaskStatus(status: Status, result: Option[ExecutionResult]) extends ExecutionEvent