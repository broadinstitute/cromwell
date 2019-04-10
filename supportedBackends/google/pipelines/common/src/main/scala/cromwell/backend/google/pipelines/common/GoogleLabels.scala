package cromwell.backend.google.pipelines.common

import cats.data.Validated.Valid
import common.validation.ErrorOr.ErrorOr
import cats.syntax.validated._
import cats.syntax.apply._

final case class GoogleLabel(key: String, value: String)

object GoogleLabels {

  val MaxLabelLength = 63
  val GoogleLabelRegexPattern = "[a-z]([-a-z0-9]*[a-z0-9])?"
  val GoogleLabelRegex = GoogleLabelRegexPattern.r

  // This function is used to coerce a string into one that meets the requirements for a label submission to Google Pipelines API.
  // See 'labels' in https://cloud.google.com/genomics/reference/rpc/google.genomics.v1alpha2#google.genomics.v1alpha2.RunPipelineArgs
  def safeGoogleName(mainText: String, emptyAllowed: Boolean = false): String = {

    validateLabelRegex(mainText) match {
      case Valid(labelText) => labelText
      case invalid @ _ if mainText.equals("") && emptyAllowed => mainText
      case invalid @ _ =>
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

  def validateLabelRegex(s: String): ErrorOr[String] = {
    (GoogleLabelRegex.pattern.matcher(s).matches, s.length <= MaxLabelLength) match {
      case (true, true) => s.validNel
      case (false, false) => s"Invalid label field: `$s` did not match regex $GoogleLabelRegexPattern and it is ${s.length} characters. The maximum is $MaxLabelLength.".invalidNel
      case (false, _) => s"Invalid label field: `$s` did not match the regex $GoogleLabelRegexPattern.".invalidNel
      case (_, false) => s"Invalid label field: `$s` is ${s.length} characters. The maximum is $MaxLabelLength.".invalidNel
    }
  }


  def safeLabels(values: (String, String)*): Seq[GoogleLabel] = {
    def safeGoogleLabel(kvp: (String, String)): GoogleLabel = {
      GoogleLabel(safeGoogleName(kvp._1), safeGoogleName(kvp._2, emptyAllowed = true))
    }
    values.map(safeGoogleLabel)
  }

  def validateLabel(key: String, value: String): ErrorOr[GoogleLabel] = {
    (validateLabelRegex(key), validateLabelRegex(value)).mapN { (validKey, validValue) => GoogleLabel(validKey, validValue)  }
  }

}
