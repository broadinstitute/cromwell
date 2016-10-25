package wdl4s.util

import java.util.regex.Pattern

import scala.annotation.tailrec

object StringUtil {
  val Ws = Pattern.compile("[\\ \\t]+")

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
    val parts = trimmed.split("[\\r\\n]+")
    val indent = parts.map(leadingWhitespaceCount).min
    parts.map(_.substring(indent)).mkString("\n")
  }

  private def leadingWhitespaceCount(s: String): Int = {
    val matcher = Ws.matcher(s)
    if (matcher.lookingAt) matcher.end else 0
  }

  def stripAll(s: String, startChars: String, endChars: String): String = {
    /* https://stackoverflow.com/questions/17995260/trimming-strings-in-scala */
    @tailrec
    def start(n: Int): String = {
      if (n == s.length) ""
      else if (startChars.indexOf(s.charAt(n)) < 0) end(n, s.length)
      else start(1 + n)
    }

    @tailrec
    def end(a: Int, n: Int): String = {
      if (n <= a) s.substring(a, n)
      else if (endChars.indexOf(s.charAt(n - 1)) < 0) s.substring(a, n)
      else end(a, n - 1)
    }

    start(0)
  }
}
