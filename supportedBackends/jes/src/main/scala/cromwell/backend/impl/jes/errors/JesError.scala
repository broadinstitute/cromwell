package cromwell.backend.impl.jes.errors

import java.nio.file.Path

import cromwell.backend.impl.jes.RunStatus

import scala.util.Try

object JesError {
  /**
    * List of known and mapped JES errors
    */
  private val KnownErrors = List(
    FailedToDelocalize
  )

  /**
    * The error message looks like "10: Something bad happened"
    * Extract 10 as an inner error code and then the rest of the message
    */
  private def splitInnerCodeAndMessage(errorMessage: String): Option[(Int, String)] = Try {
    val sep = errorMessage.indexOf(':')
    val innerCode = errorMessage.substring(0, sep).toInt
    val message = errorMessage.substring(sep + 1, errorMessage.length - 1).trim
    (innerCode, message)
  }.toOption
  
  def fromFailedStatus(failedStatus: RunStatus.Failed, jobTag: String, stderrPath: Option[Path]): Option[JesKnownJobFailure] = {
    def lookupError(innerCode: Int, message: String) = {
      val search = KnownErrors.toStream map { _.toFailureOption(failedStatus.errorCode, innerCode, message, jobTag, stderrPath) } find { _.isDefined }
      search.flatten
    }
    
    failedStatus.errorMessage flatMap splitInnerCodeAndMessage flatMap Function.tupled(lookupError)
  }
}

sealed abstract class JesError(val outerCode: Int, val innerCode: Int, val messageStart: String) {
  def toFailureOption(errorOuterCode: Int, errorInnerCode: Int, errorMessage: String, jobTag: String, stderrPath: Option[Path]): Option[JesKnownJobFailure] = {
    if (
        errorOuterCode == outerCode &&
        errorInnerCode == innerCode &&
        errorMessage.startsWith(messageStart)
    ) Option(toJobFailure(errorMessage, jobTag, stderrPath))
    else None
  }
  
  protected def toJobFailure(message: String, jobTag: String, stderrPath: Option[Path]): JesKnownJobFailure
}

private case object FailedToDelocalize extends JesError(5, 10, "Failed to delocalize files") {
  def toJobFailure(message: String, jobTag: String, stderrPath: Option[Path]) = FailedToDelocalizeFailure(message, jobTag, stderrPath)
}