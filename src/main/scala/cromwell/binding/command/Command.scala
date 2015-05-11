package cromwell.binding.command

import java.util.regex.Pattern

import cromwell.binding.WdlValue
import cromwell.binding.types.WdlType

import scala.annotation.tailrec

/**
 * Represents the `command` section of a `task` definition in a WDL file.
 *
 * A command is a sequence of CommandPart objects which can either be String
 * literals or parameters that the user needs to provide.  For example, consider
 * the following command:
 *
 * grep '${pattern}' ${file in}
 *
 * The Seq[CommandPart] would look as follows:
 *
 * <ol>
 * <li>StringCommandPart: <pre>grep '</pre></li>
 * <li>ParameterCommandPart: pattern (type=string)</li>
 * <li>StringCommandPart: <pre>' </pre></li>
 * <li>ParameterCommandPart: in (type=file)</li>
 * </ol>
 *
 * The result of the `inputs` method would be
 *
 * <ol>
 * <li>(<pre>pattern</pre>, <pre>string</pre>)</li>
 * <li>(<pre>in</pre>, <pre>file</pre>)</li>
 * </ol>
 *
 * A command line can be "instantiated" via the instantiate() method by providing a
 * value for all of its inputs.  The result is a string representation of the command
 *
 * @param parts The CommandParts that represent this command line
 */
case class Command(parts: Seq[CommandPart]) {
  val ws = Pattern.compile("[\\ \\t]+")

  def inputs: Map[String, WdlType] = parts.collect({ case p: ParameterCommandPart => (p.name, p.wdlType) }).toMap

  /**
   * Given a map of task-local parameter names and WdlValues,
   * create a command String
   *
   * @param parameters Parameter values
   * @return String instantiation of the command
   */
  def instantiate(parameters: Map[String, WdlValue]): String = normalize(parts.map { part => part.instantiate(parameters) }.mkString(""))

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
  private def normalize(s: String): String = {
    val trimmed = stripAll(s, "\r\n", "\r\n \t")
    val parts = trimmed.split("[\\r\\n]+")
    val indent = parts.map(leadingWhitespaceCount).min
    parts.map(_.substring(indent)).mkString("")
  }

  private def leadingWhitespaceCount(s: String): Int = {
    val matcher = ws.matcher(s)
    if (matcher.lookingAt) matcher.end else 0
  }

  private def stripAll(s: String, startChars: String, endChars: String): String = {
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

  override def toString: String = s"[Command: ${normalize(parts.map(y => y.toString).mkString(""))}]"
}
