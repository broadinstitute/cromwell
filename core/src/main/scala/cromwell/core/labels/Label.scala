package cromwell.core.labels

import lenthall.validation.ErrorOr.ErrorOr
import cats.data.Validated._
import cats.syntax.cartesian._
import cats.syntax.validated._

import scala.util.matching.Regex

sealed abstract case class Label(key: String, value: String)

object Label {

  // Yes, 63. Not a typo for 64.
  // See 'labels' in https://cloud.google.com/genomics/reference/rpc/google.genomics.v1alpha2#google.genomics.v1alpha2.RunPipelineArgs
  private val MaxLabelLength = 63
  val GoogleLabelRegexPattern = "[a-z]([-a-z0-9]*[a-z0-9])?"
  val LabelKeyRegex = GoogleLabelRegexPattern
  val LabelValueRegex = "([a-z0-9]*[-a-z0-9]*[a-z0-9])?"

  val LabelExpectationsMessage =
    s"A Label key must match the pattern `$LabelKeyRegex` and a label value must match the pattern `$LabelValueRegex`."

  def validateLabelRegex(s: String, regexAllowed: Regex): ErrorOr[String] = {
    (regexAllowed.pattern.matcher(s).matches, s.length <= MaxLabelLength) match {
      case (true, true) => s.validNel
      case (false, false) => s"Invalid label: `$s` did not match regex $regexAllowed and it is ${s.length} characters. The maximum is $MaxLabelLength.".invalidNel
      case (false, _) => s"Invalid label: `$s` did not match the regex $regexAllowed.".invalidNel
      case (_, false) => s"Invalid label: `$s` is ${s.length} characters. The maximum is $MaxLabelLength.".invalidNel
    }
  }

  def validateLabelKey(key: String): ErrorOr[String] = validateLabelRegex(key, LabelKeyRegex.r)

  def validateLabelValue(key: String): ErrorOr[String] = validateLabelRegex(key, LabelValueRegex.r)

  def validateLabel(key: String, value: String): ErrorOr[Label] = {
    val validatedKey = validateLabelKey(key)
    val validatedValue = validateLabelValue(value)

    (validatedKey |@| validatedValue) map { case (k, v) => new Label(k, v) {} }
  }

  /**
    * Change to meet the constraint:
    *  - To match the regex LabelRegexPattern
    *  - To be between 1 and MaxLabelLength characters total
    */
  def safeGoogleName(mainText: String, emptyAllowed: Boolean = false): String = {

    validateLabelRegex(mainText, GoogleLabelRegexPattern.r) match {
      case Valid(labelText) => labelText
      case invalid if mainText.equals("") && emptyAllowed => mainText
      case invalid =>
        def appendSafe(current: String, nextChar: Char): String = {
          nextChar match {
            case c if c.isLetterOrDigit || c == '-' => current + c.toLower
            case _ => current + '-'
          }
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
          validStartAndEnd.substring(0, subSectionLength) + middleSeparator + validStartAndEnd.substring(length - subSectionLength, length)
        } else {
          validStartAndEnd
        }
    }
  }

  def apply(key: String, value: String) = {
    new Label(key, value) {}
  }
}
