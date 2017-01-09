package cromwell.backend.impl.jes.labels

case class Label private[labels](key: String, value: String)

object Label {

  // Yes, 63. Not a typo for 64.
  // See 'labels' in https://cloud.google.com/genomics/reference/rpc/google.genomics.v1alpha2#google.genomics.v1alpha2.RunPipelineArgs
  private val MAX_LABEL_LENGTH = 63

  def validate(s: String) = "[a-z]([-a-z0-9]*[a-z0-9])?".r.pattern.matcher(s).matches && s.length <= MAX_LABEL_LENGTH

  /**
    * Change to meet the constraint:
    *  - Must match the regex [a-z]([-a-z0-9]*[a-z0-9])?
    *  - Must be between 1 and MAX_LABEL_LENGTH characters total
    */
  def safeName(mainText: String) = {

    if (validate(mainText)) {
      mainText
    } else {
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
      val tooLong = length > MAX_LABEL_LENGTH

      if (tooLong) {
        val middleSeparator = "---"
        val subSectionLength = (MAX_LABEL_LENGTH - middleSeparator.length) / 2
        validStartAndEnd.substring(0, subSectionLength) + middleSeparator + validStartAndEnd.substring(length - subSectionLength, length)
      } else {
        validStartAndEnd
      }
    }
  }

  def safeLabel(key: String, value: String): Label = {
    Label(safeName(key), safeName(value))
  }
}
