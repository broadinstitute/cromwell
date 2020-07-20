package cromwell.backend.async

import common.exception.ThrowableAggregation
import common.validation.Validation.GreaterEqualRefined
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

final case class RetryWithMoreMemory(jobTag: String, stderrPath: Option[Path]) extends KnownJobFailureException {
  override def getMessage = s"stderr for job $jobTag contained one of the `memory-retry` error-keys specified in the config. " +
    "Job might have run out of memory."
}

final case class MemoryMultiplierLessThanOne(jobTag: String,
                                             stderrPath: Option[Path],
                                             currentMultiplier: GreaterEqualRefined,
                                             memoryRetryFactor: GreaterEqualRefined) extends KnownJobFailureException {
  override def getMessage = s"The result of multiplying current memory multiplier $currentMultiplier with memory retry factor $memoryRetryFactor was not positive."
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
