package cromwell.backend.impl.jes.errors

import java.nio.file.Path

import cromwell.backend.async.KnownJobFailureException
import cromwell.core.path.PathImplicits._

sealed trait JesKnownJobFailure extends KnownJobFailureException

case class FailedToDelocalizeFailure(message: String, jobTag: String, stderrPath: Option[Path]) extends JesKnownJobFailure {
  lazy val stderrMessage = stderrPath map { p =>
    s"3) Look into the stderr (${p.toRealString}) file for evidence that some of the output files the command is expected to create were not."
  } getOrElse ""
  
  lazy val missingFilesMessage = if (message.contains("No URLs matched")) {
    s"""It appears that some of the expected output files for task $jobTag did not exist when the command exited.
       |A few things to try
       |1) Check that the output section in your WDL is correct. Remember that all output files declared in a task must exist when the command exits.
       |2) Check that the return code is available and is valid wrt your command expected exit code
       |$stderrMessage
     """.stripMargin
  } else ""
  
  override def getMessage = missingFilesMessage + message
}
