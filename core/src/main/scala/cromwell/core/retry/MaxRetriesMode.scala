package cromwell.core.retry

import scala.util.{Failure, Success, Try}

sealed trait MaxRetriesMode
case object AllErrors extends MaxRetriesMode
case object KnownErrors extends MaxRetriesMode

object MaxRetriesMode {
  val DefaultMode = AllErrors
  private val AllModes = Seq(AllErrors, KnownErrors)

  def tryParse(mode: String): Try[MaxRetriesMode] =
    AllModes find { _.toString.equalsIgnoreCase(mode) } map { Success(_) } getOrElse Failure(
      new Exception(s"Invalid max retries mode: '$mode', supported modes are: ${AllModes.mkString("'", "', '", "'")}")
    )
}
