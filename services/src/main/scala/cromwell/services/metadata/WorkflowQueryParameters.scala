package cromwell.services.metadata

import java.time.OffsetDateTime

import cats.data.Validated._
import cats.syntax.cartesian._
import cats.syntax.validated._
import cromwell.core.WorkflowId
import cromwell.services.metadata.WorkflowQueryKey._
import lenthall.validation.ErrorOr._


case class WorkflowQueryParameters private(statuses: Set[String],
                                           names: Set[String],
                                           ids: Set[WorkflowId],
                                           startDate: Option[OffsetDateTime],
                                           endDate: Option[OffsetDateTime],
                                           page: Option[Int],
                                           pageSize: Option[Int])

object WorkflowQueryParameters {

  private def validateStartBeforeEnd(start: Option[OffsetDateTime], end: Option[OffsetDateTime]): ErrorOr[Unit] = {
    // Invert the notion of success/failure here to only "successfully" generate an error message if
    // both start and end dates have been specified and start is after end.
    val startAfterEndError = for {
      s <- start
      e <- end
      if s.isAfter(e)
    } yield s"Specified start date is after specified end date: start: $s, end: $e"

    // If the Option is defined this represents a failure, if it's empty this is a success.
    startAfterEndError map { _.invalidNel } getOrElse ().validNel
  }

  private def validateOnlyRecognizedKeys(rawParameters: Seq[(String, String)]): ErrorOr[Unit] = {
    // Create a map of keys by canonical capitalization (capitalized first letter, lowercase everything else).
    // The values are the keys capitalized as actually given to the API, which is what will be used in any
    // error messages.
    val keysByCanonicalCapitalization: Map[String, Set[String]] =
      rawParameters map { _._1 } groupBy { _.toLowerCase.capitalize } mapValues { _.toSet }

    keysByCanonicalCapitalization.keys.toSet -- WorkflowQueryKey.ValidKeys match {
      case set if set.nonEmpty =>
        val unrecognized = set flatMap keysByCanonicalCapitalization
        ("Unrecognized query keys: " + unrecognized.mkString(", ")).invalidNel
      case _ => ().validNel
    }
  }

  /**
   * Run the validation logic over the specified raw parameters, creating a `WorkflowQueryParameters` if all
   * validation succeeds, otherwise accumulate all validation messages within the `ValidationNel`.
   */
  private [metadata] def runValidation(rawParameters: Seq[(String, String)]): ErrorOr[WorkflowQueryParameters] = {

    val onlyRecognizedKeys = validateOnlyRecognizedKeys(rawParameters)

    val valuesByCanonicalCapitalization: Map[String, Seq[(String, String)]] =
      rawParameters groupBy { case (key, _) => key.toLowerCase.capitalize }

    val Seq(startDate, endDate) = Seq(StartDate, EndDate) map { _.validate(valuesByCanonicalCapitalization) }
    val Seq(statuses, names, ids) = Seq(Status, Name, WorkflowQueryKey.Id) map {
      _.validate(valuesByCanonicalCapitalization)
    }

    val Seq(page, pageSize) = Seq(Page, PageSize) map { _.validate(valuesByCanonicalCapitalization) }

    // Only validate start before end if both of the individual date parsing validations have already succeeded.
    val startBeforeEnd = (startDate, endDate) match {
      case (Valid(s), Valid(e)) => validateStartBeforeEnd(s, e)
      case _ => ().validNel[String]
    }

    (onlyRecognizedKeys |@| startBeforeEnd |@| statuses |@| names |@| ids |@| startDate |@| endDate |@| page |@| pageSize) map {
      case (_, _, status, name, uuid, start, end, _page, _pageSize) =>
        val workflowId = uuid map WorkflowId.fromString
        WorkflowQueryParameters(status.toSet, name.toSet, workflowId.toSet, start, end, _page, _pageSize)
    }
  }

  def apply(rawParameters: Seq[(String, String)]): WorkflowQueryParameters = {
    runValidation(rawParameters) match {
      case Valid(queryParameters) => queryParameters
      case Invalid(x) => throw new IllegalArgumentException(x.toList.mkString("\n"))
    }
  }
}
