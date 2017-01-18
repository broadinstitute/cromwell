package cromwell.core.labels

import lenthall.validation.ErrorOr.ErrorOr
import cats.data.Validated._
import cats.syntax.cartesian._
import cats.syntax.validated._

sealed abstract case class Label(key: String, value: String)

object Label {

  // Yes, 63. Not a typo for 64.
  // See 'labels' in https://cloud.google.com/genomics/reference/rpc/google.genomics.v1alpha2#google.genomics.v1alpha2.RunPipelineArgs
  private val MaxLabelLength = 63
  val LabelRegexPattern = "[a-z]([-a-z0-9]*[a-z0-9])?"

  def validateName(s: String): ErrorOr[String] = {
    if (LabelRegexPattern.r.pattern.matcher(s).matches) {
      if (s.length <= MaxLabelLength) s.validNel else s"Invalid label: $s was ${s.length} characters. The maximum is $MaxLabelLength".invalidNel
    } else {
      s"Invalid label: $s did not match the regex $LabelRegexPattern".invalidNel
    }
  }

  def validateLabel(key: String, value: String): ErrorOr[Label] = {
    val validatedKey = validateName(key)
    val validatedValue = validateName(value)

    (validatedKey |@| validatedValue) map { case (k, v) => new Label(k, v) {} }
  }

  /**
    * Change to meet the constraint:
    *  - To match the regex LabelRegexPattern
    *  - To be between 1 and MaxLabelLength characters total
    */
  def safeName(mainText: String): String = {

    validateName(mainText) match {
      case Valid(labelText) => labelText
      case _ =>
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

  def safeLabel(key: String, value: String): Label = {
    new Label(safeName(key), safeName(value)) {}
  }
}
