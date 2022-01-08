package cromwell.backend.async

import akka.event.LoggingAdapter
import common.exception.ThrowableAggregation
import cromwell.core.path.Path
import wom.expression.{NoIoFunctionSet, WomExpression}

abstract class KnownJobFailureException extends Exception {
  def stderrPath: Option[Path]
}

final case class WrongReturnCode(jobTag: String, returnCode: Int, stderrPath: Option[Path]) extends KnownJobFailureException {
  override def getMessage = s"Job $jobTag exited with return code $returnCode which has not been declared as a valid return code. See 'continueOnReturnCode' runtime attribute for more details."
}

final case class ReturnCodeIsNotAnInt(jobTag: String, returnCode: String, stderrPath: Option[Path]) extends KnownJobFailureException {
  override def getMessage = {
    if (returnCode.isEmpty)
      s"The return code file for job $jobTag was empty."
    else
      s"Job $jobTag exited with return code $returnCode which couldn't be converted to an Integer."
  }
}

final case class StderrNonEmpty(jobTag: String, stderrLength: Long, stderrPath: Option[Path]) extends KnownJobFailureException {
  override def getMessage = s"stderr for job $jobTag has length $stderrLength and 'failOnStderr' runtime attribute was true."
}

final case class RetryWithMoreMemory(jobTag: String,
                                     stderrPath: Option[Path],
                                     memoryRetryErrorKeys: Option[List[String]],
                                     logger: LoggingAdapter) extends KnownJobFailureException {
  val errorKeysAsString = memoryRetryErrorKeys match {
    case None =>
      // this should not occur at this point as one would reach this error class only if Cromwell found one of the
      // `memory-retry-error-keys` in `stderr` of the task, which is only checked if the `memory-retry-error-keys`
      // are instantiated in Cromwell config
      logger.error(s"Programmer error: found one of the `system.memory-retry-error-keys` in the `stderr` of task but " +
        s"didn't find the error keys while generating the exception!")
      ""
    case Some(keys) => keys.mkString(": [", ",", "]")
  }
  override def getMessage = s"stderr for job `$jobTag` contained one of the `memory-retry-error-keys${errorKeysAsString}` specified in " +
    s"the Cromwell config. Job might have run out of memory."
}


object RuntimeAttributeValidationFailure {
  def apply(jobTag: String,
            runtimeAttributeName: String,
            runtimeAttributeValue: Option[WomExpression]): RuntimeAttributeValidationFailure = RuntimeAttributeValidationFailure(jobTag, runtimeAttributeName, runtimeAttributeValue, None)
}

final case class RuntimeAttributeValidationFailure private (jobTag: String,
                                                            runtimeAttributeName: String,
                                                            runtimeAttributeValue: Option[WomExpression],
                                                            stderrPath: Option[Path]) extends KnownJobFailureException {
  override def getMessage = s"Task $jobTag has an invalid runtime attribute $runtimeAttributeName = ${runtimeAttributeValue map { _.evaluateValue(Map.empty, NoIoFunctionSet)} getOrElse "!! NOT FOUND !!"}"
}

final case class RuntimeAttributeValidationFailures(throwables: List[RuntimeAttributeValidationFailure]) extends KnownJobFailureException with ThrowableAggregation {
  override def exceptionContext = "Runtime validation failed"
  override val stderrPath: Option[Path] = None
}

final case class JobAlreadyFailedInJobStore(jobTag: String, originalErrorMessage: String) extends KnownJobFailureException {
  override def stderrPath: Option[Path] = None
  override def getMessage = originalErrorMessage
}
