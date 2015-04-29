package cromwell.binding

import cromwell.parser.WdlParser._

import scala.collection.JavaConverters._

case class WdlSyntaxErrorFormatter(sourceCode: String) extends SyntaxErrorFormatter {
    def lines: Seq[String] = sourceCode.split("\n")

    private def pointToSource(line: Int, col: Int): String = s"${lines(line-1)}\n${" " * (col-1)}^"

    def unexpectedEof(method: String, expected: java.util.List[TerminalIdentifier], nt_rules: java.util.List[String]): String = "ERROR: Unexpected end of file"

    def excessTokens(method: String, terminal: Terminal): String = {
        s"""ERROR: Finished parsing without consuming all tokens.
           |
           |${pointToSource(terminal.getLine, terminal.getColumn)}
         """.stripMargin
    }

    def unexpectedSymbol(method: String, actual: Terminal, expected: java.util.List[TerminalIdentifier], rule: String): String = {
        val expectedTokens = expected.asScala.map(_.string).mkString(", ")
        s"""ERROR: Unexpected symbol (line ${actual.getLine}, col ${actual.getColumn}) when parsing '$method'.
           |
           |Expected $expectedTokens, got ${actual.toPrettyString}.
           |
           |${pointToSource(actual.getLine, actual.getColumn)}
           |
           |$rule
         """.stripMargin
    }

    def noMoreTokens(method: String, expecting: TerminalIdentifier, last: Terminal): String = {
        s"""ERROR: No more tokens.  Expecting ${expecting.string}
           |
           |${pointToSource(last.getLine, last.getColumn)}
         """.stripMargin
    }

    def invalidTerminal(method: String, invalid: Terminal): String = {
        s"""ERROR: Invalid symbol ID: ${invalid.getId} (${invalid.getTerminalStr})
           |
           |${pointToSource(invalid.getLine, invalid.getColumn)}
         """.stripMargin
    }

    def tooManyWorkflows(workflowAsts: java.util.Set[Ast]): String = {
        val otherWorkflows = workflowAsts.asScala.map({ast =>
            val name: Terminal = ast.getAttribute("name").asInstanceOf[Terminal]
            s"""Prior workflow definition (line ${name.getLine} col ${name.getColumn}):
               |
               |${pointToSource(name.getLine, name.getColumn)}
             """.stripMargin
        }).mkString("\n")

        s"""ERROR: Only one workflow definition allowed, found ${workflowAsts.size} workflows:
           |
           |$otherWorkflows
         """.stripMargin
    }

    def callReferencesBadTaskName(callAst: Ast, taskName: String): String = {
        val callTask: Terminal = callAst.getAttribute("task").asInstanceOf[Terminal]
        s"""ERROR: Call references a task ($taskName) that doesn't exist (line ${callTask.getLine}, col ${callTask.getColumn})
           |
           |${pointToSource(callTask.getLine, callTask.getColumn)}
         """.stripMargin
    }

    def callReferencesBadTaskInput(callInputAst: Ast, taskAst: Ast): String = {
        val callParameter: Terminal = callInputAst.getAttribute("key").asInstanceOf[Terminal]
        val taskName: Terminal = taskAst.getAttribute("name").asInstanceOf[Terminal]
        s"""ERROR: Call references an input on task '${taskName.getSourceString}' that doesn't exist (line ${callParameter.getLine}, col ${callParameter.getColumn})
           |
           |${pointToSource(callParameter.getLine, callParameter.getColumn)}
           |
           |Task defined here (line ${taskName.getLine}, col ${taskName.getColumn}):
           |
           |${pointToSource(taskName.getLine, taskName.getColumn)}
         """.stripMargin
    }
}
