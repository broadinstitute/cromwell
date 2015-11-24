package cromwell.webservice

import cromwell.webservice.WorkflowQueryKey._
import org.joda.time.DateTime
import org.joda.time.format.ISODateTimeFormat

import scala.language.postfixOps
import scalaz.Scalaz._
import scalaz.{Name => _, _}


case class WorkflowQueryParameters private(statuses: Set[String],
                                           names: Set[String],
                                           startDate: Option[DateTime],
                                           endDate: Option[DateTime])

object WorkflowQueryParameters {

  private lazy val Formatter = ISODateTimeFormat.dateHourMinuteSecondMillis()

  private def validateStartBeforeEnd(start: Option[DateTime], end: Option[DateTime]): ValidationNel[String, Unit] = {
    // Invert the notion of success/failure here to only "successfully" generate an error message if
    // both start and end dates have been specified and start is after end.
    val startAfterEndError = for {
      s <- start
      e <- end
      if s.isAfter(e)
    } yield s"Specified start date is after specified end date: start: ${s.toString(Formatter)}, end: ${e.toString(Formatter)}"

    // If the Option is defined this represents a failure, if it's empty this is a success.
    startAfterEndError map { _.failureNel } getOrElse ().successNel
  }

  private def validateOnlyRecognizedKeys(rawParameters: Seq[(String, String)]): ValidationNel[String, Unit] = {
    // Create a map of keys by canonical capitalization (capitalized first letter, lowercase everything else).
    // The values are the keys capitalized as actually given to the API, which is what will be used in any
    // error messages.
    val keysByCanonicalCapitalization: Map[String, Set[String]] =
      rawParameters map { _._1 } groupBy { _.toLowerCase.capitalize } mapValues { _.toSet }

    keysByCanonicalCapitalization.keys.toSet -- WorkflowQueryKey.ValidKeys match {
      case set if set.nonEmpty =>
        val unrecognized = set.flatMap { k => keysByCanonicalCapitalization.get(k).get }
        ("Unrecognized query keys: " + unrecognized.mkString(", ")).failureNel
      case _ => ().successNel
    }
  }

  /**
   * Run the validation logic over the specified raw parameters, creating a `WorkflowQueryParameters` if all
   * validation succeeds, otherwise accumulate all validation messages within the `ValidationNel`.
   */
  private [webservice] def runValidation(rawParameters: Seq[(String, String)]): ValidationNel[String, WorkflowQueryParameters] = {

    val onlyRecognizedKeys = validateOnlyRecognizedKeys(rawParameters)

    val valuesByCanonicalCapitalization: Map[String, Seq[(String, String)]] =
      rawParameters groupBy { case (key, _) => key.toLowerCase.capitalize }

    val Seq(startDate, endDate) = Seq(StartDate, EndDate) map { _.validate(valuesByCanonicalCapitalization) }
    val Seq(statuses, names) = Seq(Status, Name) map { _.validate(valuesByCanonicalCapitalization) }

    // Only validate start before end if both of the individual date parsing validations have already succeeded.
    val startBeforeEnd = (startDate, endDate) match {
      case (Success(s), Success(e)) => validateStartBeforeEnd(s, e)
      case _ => ().successNel
    }

    (onlyRecognizedKeys |@| startBeforeEnd |@| statuses |@| names |@| startDate |@| endDate) {
      case (_, _, status, name, start, end) => WorkflowQueryParameters(status.toSet, name.toSet, start, end)
    }
  }

  def apply(rawParameters: Seq[(String, String)]): WorkflowQueryParameters = {
    runValidation(rawParameters) match {
      case Success(queryParameters) => queryParameters
      case Failure(x) => throw new IllegalArgumentException(x.list.mkString("\n"))
    }
  }
}
