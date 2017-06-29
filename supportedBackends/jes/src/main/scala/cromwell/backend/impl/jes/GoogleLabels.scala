package cromwell.backend.impl.jes

import cats.data.Validated.Valid
import cromwell.core.labels.{Label, Labels}

object GoogleLabels {

  private val MaxLabelLength = Label.MaxLabelLength

  private val GoogleLabelsRegexPattern = Label.LabelKeyRegex

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
      Label(Label.safeGoogleName(key), Label.safeGoogleName(value, emptyAllowed = true))
    }

    val kvps: Seq[(String, String)] = values.toSeq
    Labels((kvps map { case(k, v) => safeGoogleLabel(k, v) } ).to[Vector])
  }
}
