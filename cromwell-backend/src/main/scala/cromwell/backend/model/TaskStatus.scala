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
  * Defines a task status with its related blob data.
  */
case class TaskStatus[T](status: Status, blob: T)
