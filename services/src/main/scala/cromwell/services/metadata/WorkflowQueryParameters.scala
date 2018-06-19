package cromwell.services.metadata

import java.time.OffsetDateTime

import cats.data.Validated._
import cats.syntax.apply._
import cats.syntax.validated._
import cromwell.core.WorkflowId
import cromwell.core.labels.Label
import cromwell.services.metadata.WorkflowQueryKey._
import common.validation.ErrorOr._

case class WorkflowQueryParameters private(statuses: Set[String],
                                           names: Set[String],
                                           ids: Set[WorkflowId],
                                           labelsAnd: Set[Label],
                                           labelsOr: Set[Label],
                                           submissionTime: Option[OffsetDateTime],
                                           startDate: Option[OffsetDateTime],
                                           endDate: Option[OffsetDateTime],
                                           page: Option[Int],
                                           pageSize: Option[Int],
                                           additionalQueryResultFields: Set[String],
                                           includeSubworkflows: Boolean)

object WorkflowQueryParameters {

  private def validateDate1BeforeDate2(date1: Option[OffsetDateTime], date2: Option[OffsetDateTime], date1Type: String, date2Type: String): ErrorOr[Unit] = {
    // Invert the notion of success/failure here to only "successfully" generate an error message if
    // both date1 and date2 dates have been specified and date1 is after date2.
    val date1AfterDate2Error = for {
      s <- date1
      e <- date2
      if s.isAfter(e)
    } yield s"Specified $date1Type date is after specified $date2Type date: $date1Type: $s, $date2Type: $e"

    // If the Option is defined this represents a failure, if it's empty this is a success.
    date1AfterDate2Error map { _.invalidNel } getOrElse ().validNel
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
    val submissionTimeValidation = SubmissionTime.validate(valuesByCanonicalCapitalization)
    val statusesValidation: ErrorOr[Set[String]] = Status.validate(valuesByCanonicalCapitalization).map(_.toSet)
    val namesValidation: ErrorOr[Set[String]] = Name.validate(valuesByCanonicalCapitalization).map(_.toSet)
    val workflowIdsValidation: ErrorOr[Set[WorkflowId]] = WorkflowQueryKey.Id.validate(valuesByCanonicalCapitalization).map(ids => (ids map WorkflowId.fromString).toSet)
    val labelsAndValidation: ErrorOr[Set[Label]] = WorkflowQueryKey.LabelAndKeyValue.validate(valuesByCanonicalCapitalization).map(_.toSet)
    val labelsOrValidation: ErrorOr[Set[Label]] = WorkflowQueryKey.LabelOrKeyValue.validate(valuesByCanonicalCapitalization).map(_.toSet)
    val pageValidation = Page.validate(valuesByCanonicalCapitalization)
    val pageSizeValidation = PageSize.validate(valuesByCanonicalCapitalization)
    val additionalQueryResultFieldsValidation: ErrorOr[Set[String]] = AdditionalQueryResultFields.validate(valuesByCanonicalCapitalization).map(_.toSet)
    val includeSubworkflowsValidation = IncludeSubworkflows.validate(valuesByCanonicalCapitalization)

    // Only validate start before end if both of the individual date parsing validations have already succeeded.
    val startBeforeEndValidation: ErrorOr[Unit] = (startDateValidation, endDateValidation) match {
      case (Valid(s), Valid(e)) => validateDate1BeforeDate2(s, e, "start", "end")
      case _ => ().validNel[String]
    }

    // Only validate submission before start if both of the individual date parsing validations have already succeeded.
    val submissionBeforeStartValidation: ErrorOr[Unit] = (submissionTimeValidation, startDateValidation) match {
      case (Valid(sub), Valid(st)) => validateDate1BeforeDate2(sub, st, "submission", "start")
      case _ => ().validNel[String]
    }

    (onlyRecognizedKeysValidation,
      startBeforeEndValidation,
      submissionBeforeStartValidation,
      statusesValidation,
      namesValidation,
      workflowIdsValidation,
      labelsAndValidation,
      labelsOrValidation,
      submissionTimeValidation,
      startDateValidation,
      endDateValidation,
      pageValidation,
      pageSizeValidation,
      additionalQueryResultFieldsValidation,
      includeSubworkflowsValidation
    ) mapN {
      (_, _, _, statuses, names, ids, labelsAnd, labelsOr, submissionTime, startDate, endDate, page, pageSize, additionalQueryResultFields, includeSubworkflows) =>
        WorkflowQueryParameters(statuses, names, ids, labelsAnd, labelsOr, submissionTime, startDate, endDate, page, pageSize, additionalQueryResultFields, includeSubworkflows)
    }
  }

  def apply(rawParameters: Seq[(String, String)]): WorkflowQueryParameters = {
    runValidation(rawParameters) match {
      case Valid(queryParameters) => queryParameters
      case Invalid(x) => throw new IllegalArgumentException(x.toList.mkString("\n"))
    }
  }
}
