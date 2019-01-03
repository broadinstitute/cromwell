package cromwell.backend.google.pipelines.common

import cats.data.Validated.Valid
import common.validation.ErrorOr.ErrorOr
import cromwell.core.labels.{Label, Labels}
import cats.syntax.validated._
import scala.util.matching.Regex

object GoogleLabels {

  val MaxLabelLength = 63
  val GoogleLabelsRegexPattern = "[a-z]([-a-z0-9]*[a-z0-9])?"

  // This function is used to coerce a string into one that meets the requirements for a label submission to JES.
  // See 'labels' in https://cloud.google.com/genomics/reference/rpc/google.genomics.v1alpha2#google.genomics.v1alpha2.RunPipelineArgs
  def safeGoogleName(mainText: String, emptyAllowed: Boolean = false): String = {

    validateLabelRegex(mainText, GoogleLabelsRegexPattern.r) match {
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

  def validateLabelRegex(s: String, regexAllowed: Regex): ErrorOr[String] = {
    (regexAllowed.pattern.matcher(s).matches, s.length <= MaxLabelLength) match {
      case (true, true) => s.validNel
      case (false, false) => s"Invalid label: `$s` did not match regex $regexAllowed and it is ${s.length} characters. The maximum is $MaxLabelLength.".invalidNel
      case (false, _) => s"Invalid label: `$s` did not match the regex $regexAllowed.".invalidNel
      case (_, false) => s"Invalid label: `$s` is ${s.length} characters. The maximum is $MaxLabelLength.".invalidNel
    }
  }


  def toLabels(values: (String, String)*): Labels = {

    def safeGoogleLabel(key: String, value: String): Label = {
      Label(safeGoogleName(key), safeGoogleName(value, emptyAllowed = true))
    }

    Labels(values.toVector map (safeGoogleLabel _ ).tupled)
  }
}
