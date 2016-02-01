package cromwell.backend.model

import wdl4s.values.WdlValue

/**
  * Represents a result of task execution.
  */
trait ExecutionResult

/**
  * Successful task execution.
  * @param outputs Outputs from task.
  */
case class SuccessfulTaskResult(outputs: Map[String, WdlValue]) extends ExecutionResult

/**
  * Failure task execution.
  * @param exception Exception generated during task execution.
  * @param rc Return code from command.
  * @param stdErr Standard error data.
  */
case class FailureTaskResult(exception: Throwable, rc: Int, stdErr: String) extends ExecutionResult

/**
  * Failure execution.
  * @param exception Exception generated during task execution.
  */
case class FailureResult(exception: Throwable) extends ExecutionResult