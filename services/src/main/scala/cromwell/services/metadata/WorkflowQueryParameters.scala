package cromwell.services.metadata

import java.time.OffsetDateTime

import cats.data.Validated._
import cats.syntax.apply._
import cats.syntax.validated._
import cromwell.core.WorkflowId
import cromwell.core.labels.Label
import cromwell.services.metadata.WorkflowQueryKey._
import lenthall.validation.ErrorOr._

case class WorkflowQueryParameters private(statuses: Set[String],
                                           names: Set[String],
                                           ids: Set[WorkflowId],
                                           labels: Set[Label],
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

    val onlyRecognizedKeysValidation = validateOnlyRecognizedKeys(rawParameters)

    val valuesByCanonicalCapitalization: Map[String, Seq[(String, String)]] =
      rawParameters groupBy { case (key, _) => key.toLowerCase.capitalize }

    val startDateValidation = StartDate.validate(valuesByCanonicalCapitalization)
    val endDateValidation = EndDate.validate(valuesByCanonicalCapitalization)
    val statusesValidation: ErrorOr[Set[String]] = Status.validate(valuesByCanonicalCapitalization).map(_.toSet)
    val namesValidation: ErrorOr[Set[String]] = Name.validate(valuesByCanonicalCapitalization).map(_.toSet)
    val workflowIdsValidation: ErrorOr[Set[WorkflowId]] = WorkflowQueryKey.Id.validate(valuesByCanonicalCapitalization).map(ids => (ids map WorkflowId.fromString).toSet)
    val labelsValidation: ErrorOr[Set[Label]] = WorkflowQueryKey.LabelKeyValue.validate(valuesByCanonicalCapitalization).map(_.toSet)
    val pageValidation = Page.validate(valuesByCanonicalCapitalization)
    val pageSizeValidation = PageSize.validate(valuesByCanonicalCapitalization)

    // Only validate start before end if both of the individual date parsing validations have already succeeded.
    val startBeforeEndValidation: ErrorOr[Unit] = (startDateValidation, endDateValidation) match {
      case (Valid(s), Valid(e)) => validateStartBeforeEnd(s, e)
      case _ => ().validNel[String]
    }

    (onlyRecognizedKeysValidation,
      startBeforeEndValidation,
      statusesValidation,
      namesValidation,
      workflowIdsValidation,
      labelsValidation,
      startDateValidation,
      endDateValidation,
      pageValidation,
      pageSizeValidation) mapN { (_, _, statuses, names, ids, labels, startDate, endDate, page, pageSize) => WorkflowQueryParameters(statuses, names, ids, labels, startDate, endDate, page, pageSize) }
  }

  def apply(rawParameters: Seq[(String, String)]): WorkflowQueryParameters = {
    runValidation(rawParameters) match {
      case Valid(queryParameters) => queryParameters
      case Invalid(x) => throw new IllegalArgumentException(x.toList.mkString("\n"))
    }
  }
}
