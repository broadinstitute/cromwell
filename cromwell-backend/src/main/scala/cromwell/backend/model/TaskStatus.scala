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
  * Defines a task intermediate status.
  */
case class TaskStatus(status: Status) extends ExecutionEvent

/**
  * Defines a task final status with the resulting data.
  */
case class TaskFinalStatus(status: Status, result: ExecutionResult) extends ExecutionEvent