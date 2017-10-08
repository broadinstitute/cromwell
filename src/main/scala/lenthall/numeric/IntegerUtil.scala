package lenthall.numeric

object IntegerUtil {

  private def ordinal(int: Int): String = {
    int match {
      case 1 => "st"
      case 2 => "nd"
      case 3 => "rd"
      case _ => "th"
    }
  }

  implicit class IntEnhanced(val value: Int) extends AnyVal {
    def toOrdinal: String = value match {
      case v if v.isBetweenInclusive(10, 20) => s"${v}th"
      case v =>
        val suffix = ordinal(v % 10)
        s"$v$suffix"
    }

    def isBetweenInclusive(min: Int, max: Int): Boolean = {
      min <= value && value <= max
    }
  }

}
