package cromwell.binding

import cromwell.parser.WdlParser._

import scala.collection.JavaConverters._

case class WdlSyntaxErrorFormatter(terminalMap: Map[Terminal, WdlSource]) extends SyntaxErrorFormatter {

  private def pointToSource(t: Terminal): String = s"${line(t)}\n${" " * (t.getColumn - 1)}^"
  private def line(t:Terminal): String = terminalMap.get(t).get.split("\n")(t.getLine - 1)

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
        |Expected $expectedTokens, got ${actual.toPrettyString}.
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

  def tooManyWorkflows(workflowAsts: java.util.List[Ast]): String = {
    val otherWorkflows = workflowAsts.asScala.map({ ast =>
      val name: Terminal = ast.getAttribute("name").asInstanceOf[Terminal]
      s"""Prior workflow definition (line ${name.getLine} col ${name.getColumn}):
          |
          |${pointToSource(name)}
       """.stripMargin
    }).mkString("\n")

    s"""ERROR: Only one workflow definition allowed, found ${workflowAsts.size} workflows:
        |
        |$otherWorkflows
     """.stripMargin
  }

  def duplicateTask(taskAsts: Seq[Ast]): String = {
    val otherTasks = taskAsts.map({ ast =>
      val name: Terminal = ast.getAttribute("name").asInstanceOf[Terminal]
      s"""Prior task definition (line ${name.getLine} col ${name.getColumn}):
       |
       |${pointToSource(name)}
       """.stripMargin
    }).mkString("\n")

    s"""ERROR: Two tasks defined with the name '${taskAsts.head.getAttribute("name").asInstanceOf[Terminal].getSourceString}':
     |
     |$otherTasks
     """.stripMargin
  }

  def callReferencesBadTaskName(callAst: Ast, taskName: String): String = {
    val callTask: Terminal = callAst.getAttribute("task").asInstanceOf[Terminal]
    s"""ERROR: Call references a task ($taskName) that doesn't exist (line ${callTask.getLine}, col ${callTask.getColumn})
        |
        |${pointToSource(callTask)}
     """.stripMargin
  }

  def callReferencesBadTaskInput(callInputAst: Ast, taskAst: Ast): String = {
    val callParameter: Terminal = callInputAst.getAttribute("key").asInstanceOf[Terminal]
    val taskName: Terminal = taskAst.getAttribute("name").asInstanceOf[Terminal]
    s"""ERROR: Call references an input on task '${taskName.getSourceString}' that doesn't exist (line ${callParameter.getLine}, col ${callParameter.getColumn})
        |
        |${pointToSource(callParameter)}
        |
        |Task defined here (line ${taskName.getLine}, col ${taskName.getColumn}):
        |
        |${pointToSource(taskName)}
     """.stripMargin
  }

  def taskHasDuplicatedInputs(taskAst: Ast, duplicatedInputAsts: Seq[Ast]): String = {
    val inputName = duplicatedInputAsts.head.getAttribute("name").asInstanceOf[Terminal].getSourceString
    val duplicatedInputs = duplicatedInputAsts.map({ ast =>
      val name: Terminal = ast.getAttribute("name").asInstanceOf[Terminal]
      s"""Input defined here (line ${name.getLine} col ${name.getColumn}):
         |
         |${pointToSource(name)}
       """.stripMargin
    }).mkString("\n")
    val taskName: Terminal = taskAst.getAttribute("name").asInstanceOf[Terminal]
    s"""ERROR: Task '${taskName.getSourceString}' has duplicated input '$inputName':
       |
       |$duplicatedInputs
       |
       |Task defined here (line ${taskName.getLine}, col ${taskName.getColumn}):
       |
       |${pointToSource(taskName)}
     """.stripMargin
  }

  def taskAndNamespaceHaveSameName(taskAst: Ast, namespace: Terminal): String = {
    val taskName = taskAst.getAttribute("name").asInstanceOf[Terminal]
    s"""ERROR: Task and namespace have the same name:
     |
     |Task defined here (line ${taskName.getLine}, col ${taskName.getColumn}):
     |
     |${pointToSource(taskName)}
     |
     |Import statement defined here (line ${namespace.getLine}, col ${namespace.getColumn}):
     |
     |${pointToSource(namespace)}
     """.stripMargin
  }

  def workflowAndNamespaceHaveSameName(workflowAst: Ast, namespace: Terminal): String = {
    val workflowName = workflowAst.getAttribute("name").asInstanceOf[Terminal]
    s"""ERROR: Task and namespace have the same name:
     |
     |Task defined here (line ${workflowName.getLine}, col ${workflowName.getColumn}):
     |
     |${pointToSource(workflowName)}
     |
     |Import statement defined here (line ${namespace.getLine}, col ${namespace.getColumn}):
     |
     |${pointToSource(namespace)}
     """.stripMargin
  }

  def undefinedMemberAccess(ast: Ast): String = {
    val rhsAst = ast.getAttribute("rhs").asInstanceOf[Terminal]
    s"""ERROR: Expression will not evaluate (line ${rhsAst.getLine}, col ${rhsAst.getColumn}):
     |
     |${pointToSource(rhsAst)}
     """.stripMargin
  }

  def memberAccessReferencesBadTaskInput(ast: Ast): String = {
    val rhsAst = ast.getAttribute("rhs").asInstanceOf[Terminal]
    s"""ERROR: Expression reference input on task that doesn't exist (line ${rhsAst.getLine}, col ${rhsAst.getColumn}):
     |
     |${pointToSource(rhsAst)}
     """.stripMargin
  }

  def arrayMustHaveOnlyOneTypeParameter(arrayDecl: Terminal): String = {
    s"""ERROR: Array type should only have one parameterized type (line ${arrayDecl.getLine}, col ${arrayDecl.getColumn}):
     |
     |${pointToSource(arrayDecl)}
     """.stripMargin
  }

  def arrayMustHaveATypeParameter(arrayDecl: Terminal): String = {
    s"""ERROR: Array type should have exactly one parameterized type (line ${arrayDecl.getLine}, col ${arrayDecl.getColumn}):
     |
     |${pointToSource(arrayDecl)}
     """.stripMargin
  }

  def postfixQualifierRequiresSeparator(quantifier: Terminal) = {
    s"""ERROR: Parameters that specify * or + must also specify sep=""
       |
       |${pointToSource(quantifier)}
     """.stripMargin
  }
}
