package wdl.draft3.transforms.parsing

import scala.collection.JavaConverters._
import wdl.draft3.parser.WdlParser._
import wom.core.WorkflowSource

case class WdlDraft3SyntaxErrorFormatter(terminalMap: Map[Terminal, WorkflowSource]) extends SyntaxErrorFormatter {

  private def pointToSource(t: Terminal): String = s"${line(t)}\n${" " * (t.getColumn - 1)}^"

  private def getTerminal(t: Terminal) = t match {
    case classicTerminal => terminalMap.get(classicTerminal)
  }

  private def line(t: Terminal): String = getTerminal(t).map(_.split("\n")(t.getLine - 1)).getOrElse(s"Cannot highlight line. It was probably in an imported file.")

  def unexpectedEof(method: String, expected: java.util.List[TerminalIdentifier], nt_rules: java.util.List[String]): String = "ERROR: Unexpected end of file"

  def excessTokens(method: String, terminal: Terminal): String = {
    s"""ERROR: Finished parsing without consuming all tokens.
       |
        |${pointToSource(terminal)}
     """.stripMargin
  }

  def unexpectedSymbol(method: String, actual: Terminal, expected: java.util.List[TerminalIdentifier], rule: String): String = {
    val expectedTokens = expected.asScala.map(_.string).mkString(", ")
    s"""ERROR: Unexpected symbol (line ${actual.getLine}, col ${actual.getColumn}) when parsing '$method'.
       |
        |Expected $expectedTokens, got "${actual.getSourceString}".
       |
        |${pointToSource(actual)}
       |
        |$rule
     """.stripMargin
  }

  def noMoreTokens(method: String, expecting: TerminalIdentifier, last: Terminal): String = {
    s"""ERROR: No more tokens.  Expecting ${expecting.string}
       |
        |${pointToSource(last)}
     """.stripMargin
  }

  def invalidTerminal(method: String, invalid: Terminal): String = {
    s"""ERROR: Invalid symbol ID: ${invalid.getId} (${invalid.getTerminalStr})
       |
        |${pointToSource(invalid)}
     """.stripMargin
  }

  // TODO: these next two methods won't be called by the parser because there are no lists in the WDL grammar that
  // cause these to be triggered.  Currently the parser is passing in 'null' for the value of 'last' and when that
  // changes, these errors can be made more helpful.

  def missingListItems(method: String, required: Int, found: Int, last: Terminal): String = {
    s"ERROR: $method requires $required items, but only found $found"
  }

  def missingTerminator(method: String, terminal: TerminalIdentifier, last: Terminal): String = {
    s"ERROR: $method requires a terminator after each element"
  }
}
