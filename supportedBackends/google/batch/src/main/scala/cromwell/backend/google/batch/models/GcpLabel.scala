package cromwell.backend.google.batch.models

import cats.data.Validated.{Invalid, Valid}
import cats.instances.list._
import cats.syntax.apply._
import cats.syntax.traverse._
import cats.syntax.validated._
import common.exception.AggregatedMessageException
import common.validation.ErrorOr.ErrorOr
import cromwell.core.{CromwellFatalExceptionMarker, WorkflowOptions}
import spray.json.{JsObject, JsString}

import scala.util.control.NoStackTrace
import scala.util.{Failure, Success, Try}

final case class GcpLabel(key: String, value: String)

object GcpLabel {

  val MaxLabelLength = 63
  val GoogleLabelKeyRegexPattern = "[a-z]([-a-z_0-9]*[a-z0-9])?" // label key must start with letter
  val GoogleLabelValueRegexPattern = "[a-z0-9]([-a-z_0-9]*[a-z0-9])?"
  val GoogleLabelKeyRegex = GoogleLabelKeyRegexPattern.r
  val GoogleLabelValueRegex = GoogleLabelValueRegexPattern.r

  // This function is used to coerce a string into one that meets the requirements for a label submission to Google Batch API.
  // https://cloud.google.com/compute/docs/labeling-resources#requirements
  def safeGoogleName(mainText: String, emptyAllowed: Boolean = false): String =
    validateLabelRegex(mainText, true) match {
      case Valid(labelText) => labelText
      case invalid @ _ if mainText.equals("") && emptyAllowed => mainText
      case invalid @ _ =>
        def appendSafe(current: String, nextChar: Char): String =
          nextChar match {
            case c if c.isLetterOrDigit || c == '-' => current + c.toLower
            case _ => current + '-'
          }

        val foldResult = mainText.toCharArray.foldLeft("")(appendSafe)

        val startsValid = foldResult.headOption.exists(_.isLetter)
        val endsValid = foldResult.lastOption.exists(_.isLetterOrDigit)

        val validStart = if (startsValid) foldResult else "x--" + foldResult
        val validStartAndEnd = if (endsValid) validStart else validStart + "--x"

        val length = validStartAndEnd.length
        val tooLong = length > MaxLabelLength

        if (tooLong) {
          val middleSeparator = "---"
          val subSectionLength = (MaxLabelLength - middleSeparator.length) / 2
          validStartAndEnd.substring(0, subSectionLength) + middleSeparator + validStartAndEnd.substring(
            length - subSectionLength,
            length
          )
        } else {
          validStartAndEnd
        }
    }

  def validateLabelRegex(s: String, isKey: Boolean): ErrorOr[String] = {
    val regexPattern = if (isKey) GoogleLabelKeyRegex.pattern else GoogleLabelValueRegex.pattern
    val field = if (isKey) "label key" else "label value"

    (regexPattern.matcher(s).matches, s.length <= MaxLabelLength) match {
      case (true, true) => s.validNel
      case (false, false) =>
        s"Invalid $field field: `$s` did not match regex '$regexPattern' and it is ${s.length} characters. The maximum is $MaxLabelLength.".invalidNel
      case (false, _) => s"Invalid $field field: `$s` did not match the regex '$regexPattern'".invalidNel
      case (_, false) =>
        s"Invalid $field field: `$s` is ${s.length} characters. The maximum is $MaxLabelLength.".invalidNel
    }
  }

  def safeLabels(values: (String, String)*): Seq[GcpLabel] = {
    def safeGoogleLabel(kvp: (String, String)): GcpLabel =
      GcpLabel(safeGoogleName(kvp._1), safeGoogleName(kvp._2, emptyAllowed = true))
    values.map(safeGoogleLabel)
  }

  def validateLabel(key: String, value: String): ErrorOr[GcpLabel] =
    (validateLabelRegex(key, true), validateLabelRegex(value, false)).mapN { (validKey, validValue) =>
      GcpLabel(validKey, validValue)
    }

  def fromWorkflowOptions(workflowOptions: WorkflowOptions): Try[Seq[GcpLabel]] = {

    def extractGoogleLabelsFromJsObject(jsObject: JsObject): Try[Seq[GcpLabel]] = {
      val asErrorOr = jsObject.fields.toList.traverse {
        case (key: String, value: JsString) => GcpLabel.validateLabel(key, value.value)
        case (key, other) =>
          s"Bad label value type for '$key'. Expected simple string but got $other".invalidNel: ErrorOr[GcpLabel]
      }

      asErrorOr match {
        case Valid(value) => Success(value)
        case Invalid(errors) =>
          Failure(
            new AggregatedMessageException("Invalid 'google_labels' in workflow options", errors.toList)
              with CromwellFatalExceptionMarker
              with NoStackTrace
          )
      }
    }

    workflowOptions.toMap.get("google_labels") match {
      case Some(obj: JsObject) => extractGoogleLabelsFromJsObject(obj)
      case Some(other) =>
        Failure(
          new Exception(
            s"Invalid 'google_labels' in workflow options. Must be a simple JSON object mapping string keys to string values. Got $other"
          ) with NoStackTrace with CromwellFatalExceptionMarker
        )
      case None => Success(Seq.empty)
    }
  }

}
