package cromwell.backend.async

import java.nio.file.Path

sealed abstract class KnownJobFailureException extends Exception {
  def stderrPath: Path
  def jobTag: String
}

final case class WrongReturnCode(jobTag: String, returnCode: Int, stderrPath: Path) extends KnownJobFailureException {
  override def getMessage = s"Job $jobTag exited with return code $returnCode which has not been declared as a valid return code. See 'continueOnReturnCode' runtime attribute for more details."
}

final case class ReturnCodeIsNotAnInt(jobTag: String, returnCode: String, stderrPath: Path) extends KnownJobFailureException {
  override def getMessage = {
    if (returnCode.isEmpty)
      s"The return code file for job $jobTag was empty."
    else
      s"Job $jobTag exited with return code $returnCode which couldn't be converted to an Integer."
  }
}

final case class StderrNonEmpty(jobTag: String, stderrLength: Long, stderrPath: Path) extends KnownJobFailureException {
  override def getMessage = s"stderr for job $jobTag has length $stderrLength and 'failOnStderr' runtime attribute was true."
}
