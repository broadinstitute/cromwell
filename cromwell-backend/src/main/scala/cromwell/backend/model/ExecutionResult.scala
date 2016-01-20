package cromwell.backend.model

import wdl4s.values.WdlValue

/**
  * Represents a result of task execution.
  */
trait ExecutionResult

/**
  * Succesful execution.
  * @param outputs Outputs from task.
  */
case class SuccesfulResult(outputs: Map[String, WdlValue]) extends ExecutionResult

/**
  * Failure execution.
  * @param exception Exception generated during task execution.
  * @param rc Return code from command.
  * @param stdErr Standard error data.
  */
case class FailureResult(exception: Throwable, rc: Option[Int], stdErr: String) extends ExecutionResult
