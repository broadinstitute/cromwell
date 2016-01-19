package cromwell.backend.model

import wdl4s.values.WdlValue

trait ExecutionResult

case class SuccesfulResult(outputs: Map[String, WdlValue]) extends ExecutionResult

case class FailureResult(exception: Throwable, rc: Option[Int], stdErr: String) extends ExecutionResult
