package cromwell.backend.google.pipelines.common.errors

import cromwell.backend.async.KnownJobFailureException
import cromwell.core.path.Path

sealed trait PipelinesApiKnownJobFailure extends KnownJobFailureException

case class FailToLocalizeFailure(messages: List[String]) extends PipelinesApiKnownJobFailure {
  override def stderrPath = None
}

case class FailedToDelocalizeFailure(message: String, jobTag: String, stderrPath: Option[Path])
  extends PipelinesApiKnownJobFailure {
  lazy val stderrMessage = stderrPath map { p =>
    s"3) Look into the stderr (${p.pathAsString}) file for evidence that some of the output files the command is expected to create were not created."
  } getOrElse ""

  lazy val missingFilesMessage = if (message.contains("No URLs matched")) {
    s"""It appears that some of the expected output files for task $jobTag did not exist when the command exited.
       |A few things to try
       |1) Check that the output section in your workflow is correct. Remember that all non optional output files declared in a task must exist when the command exits.
       |2) Check that the return code is available and is valid with respect to your command expected exit code
       |$stderrMessage
     """.stripMargin
  } else ""

  override def getMessage = missingFilesMessage + message
}
