package wdl.util

import java.util.regex.Pattern

import scala.annotation.tailrec

/** This is under a `wdl` package because it exists purely to do WDLy stuff, but it's currently being called from the
  * WOMmy TaskDefinition. That should get straightened out. */
object StringUtil {
  val Ws = Pattern.compile("[\\ \\t]+")
  val utf8mb4Regex = "[\\x{10000}-\\x{FFFFF}]"
  val utf8mb3Replacement = "\uFFFD" // This is the standard char for replacing

  /**
    * 1) Remove all leading newline chars
    * 2) Remove all trailing newline AND whitespace chars
    * 3) Remove all *leading* whitespace that's common among every line in the input string
    *
    * For example, the input string:
    *
    * "
    *   first line
    *     second line
    *   third line
    *
    * "
    *
    * Would be normalized to:
    *
    * "first line
    *   second line
    * third line"
    *
    * @param s String to process
    * @return String which has common leading whitespace removed from each line
    */
  def normalize(s: String): String = {
    val trimmed = stripAll(s, "\r\n", "\r\n \t")
    val parts = trimmed.split("\\r?\\n")
    val indent = parts.filterNot(_.trim.isEmpty).map(leadingWhitespaceCount).toList match {
      case Nil => 0
      case nonEmpty => nonEmpty.min
    }
    parts.map(_.drop(indent)).mkString("\n")
  }

  private def leadingWhitespaceCount(s: String): Int = {
    val matcher = Ws.matcher(s)
    if (matcher.lookingAt) matcher.end else 0
  }

  def stripAll(s: String, startChars: String, endChars: String): String = {
    /* https://stackoverflow.com/questions/17995260/trimming-strings-in-scala */
    @tailrec
    def start(n: Int): String =
      if (n == s.length) ""
      else if (startChars.indexOf(s.charAt(n).toInt) < 0) end(n, s.length)
      else start(1 + n)

    @tailrec
    def end(a: Int, n: Int): String =
      if (n <= a) s.substring(a, n)
      else if (endChars.indexOf(s.charAt(n - 1).toInt) < 0) s.substring(a, n)
      else end(a, n - 1)

    start(0)
  }

  /**
   * Remove all utf8mb4 exclusive characters (emoji) from the given string.
   * @param in String to clean
   * @return String with all utf8mb4 exclusive characters removed
   */
  def cleanUtf8mb4(in: String): String =
    in.replaceAll(utf8mb4Regex, utf8mb3Replacement)
}
