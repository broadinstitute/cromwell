package cromwell.backend.impl.jes

import cats.data.Validated.Valid
import cromwell.core.labels.{Label, Labels}

object GoogleLabels {

  val MaxLabelLength = Label.MaxLabelLength
  val GoogleLabelsRegexPattern = Label.LabelKeyRegex

  // This function is used to coerce a string into one that meets the requirements for a label submission to JES.
  // See 'labels' in https://cloud.google.com/genomics/reference/rpc/google.genomics.v1alpha2#google.genomics.v1alpha2.RunPipelineArgs
  def safeGoogleName(mainText: String, emptyAllowed: Boolean = false): String = {

    Label.validateLabelRegex(mainText, GoogleLabelsRegexPattern.r) match {
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

  def toLabels(values: (String, String)*): Labels = {

    def safeGoogleLabel(key: String, value: String): Label = {
      Label(safeGoogleName(key), safeGoogleName(value, emptyAllowed = true))
    }

    Labels(values.toVector map (safeGoogleLabel _ ).tupled)
  }
}
