package wdl4s

import wdl4s.types.WdlType
import wdl4s.parser.WdlParser._

import scala.collection.JavaConverters._

case class WdlSyntaxErrorFormatter(terminalMap: Map[Terminal, WdlSource]) extends SyntaxErrorFormatter {

  private def pointToSource(t: Terminal): String = s"${line(t)}\n${" " * (t.getColumn - 1)}^"
  private def line(t:Terminal): String = terminalMap(t).split("\n")(t.getLine - 1)

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
        |Expected $expectedTokens, got ${actual.getSourceString}.
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

  def callReferencesAbsentTaskInput(callInputAst: Ast, taskAst: Ast, missingInput: String, callName: String): String = {
    val callParameter: Terminal = callInputAst.getAttribute("key").asInstanceOf[Terminal]
    val taskName: Terminal = taskAst.getAttribute("name").asInstanceOf[Terminal]
    val taskNameString = taskName.getSourceString
    s"""ERROR: Call supplied an unexpected input: The '$taskNameString' task doesn't have an input called '$missingInput':
        |
        |${pointToSource(callParameter)}
        |
        |Options:
        | - Add the input '$missingInput' to the '$taskNameString' task (defined on line ${taskName.getLine}).
        | - Remove '$missingInput = ...' from $callName's inputs (on line ${callParameter.getLine}).
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
    s"""ERROR: Workflow and namespace have the same name:
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

  def twoSiblingScopesHaveTheSameName(firstScopeType: String, firstScopeName: Terminal, secondScopeType: String, secondScopeName: Terminal): String = {
    s"""ERROR: Sibling nodes have conflicting names:
       |
       |$firstScopeType defined here (line ${firstScopeName.getLine}, col ${firstScopeName.getColumn}):
       |
       |${pointToSource(firstScopeName)}
       |
       |$secondScopeType statement defined here (line ${secondScopeName.getLine}, col ${secondScopeName.getColumn}):
       |
       |${pointToSource(secondScopeName)}
     """.stripMargin
  }

  def multipleCallsAndHaveSameName(names: Seq[(String, Terminal)]): String = {
    val duplicatedCallNames = names.map { case (astType, name) =>
      s"""${astType} statement here (line ${name.getLine}, column ${name.getColumn}):
        |
        |${pointToSource(name)}
      """.stripMargin
    }

    s"""ERROR: Two or more calls or values in the workflow have the same name:
       |
       |${duplicatedCallNames.mkString("\n")}
     """.stripMargin
  }

  def multipleInputStatementsOnCall(secondInputStatement: Terminal): String = {
    s"""ERROR: Call has multiple 'input' sections defined:
       |
       |${pointToSource(secondInputStatement)}
       |
       |Instead of multiple 'input' sections, use commas to separate the values.
     """.stripMargin
  }

  def emptyInputSection(callTaskName: Terminal) = {
    s"""ERROR: empty "input" section for call '${callTaskName.getSourceString}':
       |
       |${pointToSource(callTaskName)}
     """.stripMargin
  }

  def undefinedMemberAccess(ast: Ast): String = {
    val rhsAst = ast.getAttribute("rhs").asInstanceOf[Terminal]
    s"""ERROR: Expression will not evaluate (line ${rhsAst.getLine}, col ${rhsAst.getColumn}):
     |
     |${pointToSource(rhsAst)}
     """.stripMargin
  }

  def memberAccessReferencesAbsentCallOutput(memberAccessAst: Ast, call: Call): String = {
    val rhsAst = memberAccessAst.getAttribute("rhs").asInstanceOf[Terminal]
    val memberAccess = MemberAccess(memberAccessAst)
    val taskName = call.unqualifiedName
    val goodOutputs = s" (current outputs of '$taskName': " + call.outputs.map("'" + _.unqualifiedName + "'").mkString(", ") + ")"

    s"""ERROR: Call output not found: Call '${memberAccess.lhs}' doesn't have an output '${memberAccess.rhs}' (line ${rhsAst.getLine}, col ${rhsAst.getColumn}).
     |
     |${pointToSource(rhsAst)}
     |
     |Options:
     | - Add the output '${memberAccess.rhs}' to '$taskName'.
     | - Modify the member access (on line ${rhsAst.getLine}) to use an existing output$goodOutputs.
     """.stripMargin
  }

  def badOldStyleWorkflowOutput(ast: Ast): String = {
    val rhsAst = ast.getAttribute("fqn").asInstanceOf[Terminal]

    s"""ERROR: Old style workflow output references '${rhsAst.getSourceString}' which doesn't exist (line ${rhsAst.getLine}, col ${rhsAst.getColumn}):
        |
        |${pointToSource(rhsAst)}
        |""".stripMargin
  }

  def variableIsNotAnObject(ast: Ast): String = {
    val lhsAst = ast.getAttribute("lhs").asInstanceOf[Terminal]
    s"""ERROR: Variable is not an object (line ${lhsAst.getLine}, col ${lhsAst.getColumn}):
        |
     |${pointToSource(lhsAst)}
     """.stripMargin
  }

  def pairMustHaveExactlyTwoTypeParameters(arrayDecl: Terminal): String = {
    s"""ERROR: Pair type should have exactly two parameterized types (line ${arrayDecl.getLine}, col ${arrayDecl.getColumn}):
        |
     |${pointToSource(arrayDecl)}
     """.stripMargin
  }

  def arrayMustHaveOnlyOneTypeParameter(arrayDecl: Terminal): String = {
    s"""ERROR: Array type should only have one parameterized type (line ${arrayDecl.getLine}, col ${arrayDecl.getColumn}):
     |
     |${pointToSource(arrayDecl)}
     """.stripMargin
  }

  def mapMustHaveExactlyTwoTypeParameters(mapDecl: Terminal): String = {
    s"""ERROR: Map type should have two parameterized types (line ${mapDecl.getLine}, col ${mapDecl.getColumn}):
     |
     |${pointToSource(mapDecl)}
     """.stripMargin
  }

  def arrayMustHaveATypeParameter(arrayDecl: Terminal): String = {
    s"""ERROR: Array type should have exactly one parameterized type (line ${arrayDecl.getLine}, col ${arrayDecl.getColumn}):
     |
     |${pointToSource(arrayDecl)}
     """.stripMargin
  }

  def taskOutputExpressionTypeDoesNotMatchDeclaredType(outputName: Terminal, outputType: WdlType, expressionType: WdlType) = {
    s"""ERROR: ${outputName.getSourceString} is declared as a ${outputType.toWdlString} but the expression evaluates to a ${expressionType.toWdlString}:
       |
       |${pointToSource(outputName)}
     """.stripMargin
  }

  def declarationExpressionNotCoerceableToTargetType(declName: Terminal, declType: WdlType) = {
    s"""ERROR: Value for ${declName.getSourceString} is not coerceable into a ${declType.toWdlString}:
        |
       |${pointToSource(declName)}
     """.stripMargin
  }

  def failedToDetermineTypeOfDeclaration(declName: Terminal) = {
    s"""ERROR: Could not determine type of declaration ${declName.getSourceString}:
       |
       |${pointToSource(declName)}
     """.stripMargin
  }

  def trueAndFalseAttributesAreRequired(firstAttribute: Terminal) = {
    s"""ERROR: Both 'true' and 'false' attributes must be specified if either is specified:
       |
       |${pointToSource(firstAttribute)}
     """.stripMargin
  }

  def expressionExpectedToBeString(key: Terminal) = {
    s"""ERROR: Value for this attribute is expected to be a string:
       |
       |${pointToSource(key)}
     """.stripMargin
  }

  def expectedAtMostOneSectionPerTask(section: String, taskName: Terminal) = {
    s"""ERROR: Expecting to find at most one '$section' section in the task:
       |
       |${pointToSource(taskName)}
     """.stripMargin
  }

  def expectedExactlyOneCommandSectionPerTask(taskName: Terminal) = {
    s"""ERROR: Expecting to find at most one 'command' section in the task:
       |
       |${pointToSource(taskName)}
     """.stripMargin
  }

  def commandExpressionContainsInvalidVariableReference(taskName: Terminal, variable: Terminal) = {
    s"""ERROR: Variable ${variable.getSourceString} does not reference any declaration in the task (line ${variable.getLine}, col ${variable.getColumn}):
        |
        |${pointToSource(variable)}
        |
        |Task defined here (line ${taskName.getLine}, col ${taskName.getColumn}):
        |
        |${pointToSource(taskName)}
     """.stripMargin
  }

  def declarationContainsReferenceToAbsentValue(parent: Option[Scope], variable: Terminal) = {

    val (parentName, missingType) = parent match {
      case Some(t: Task) => (s"task '${t.unqualifiedName}'" , "value")
      case Some(t: TaskCall) => (s"task '${t.task.unqualifiedName}'" , "value")
      case Some(other) => ("workflow", "value or call")
      case None => ("", "")
    }

    invalidVariableReference(variable, missingType, parentName)
  }

  def scatterCollectionContainsInvalidVariableReference(scatter: Scatter, variable: Terminal) = invalidVariableReference(variable, "value", "workflow")

  private def invalidVariableReference(variable: Terminal, missingType: String, parentName: String) =
    s"""
       |ERROR: Missing $missingType: Couldn't find $missingType with name '${variable.getSourceString}' in $parentName (line ${variable.getLine}):
       |
       |${pointToSource(variable)}
       |""".stripMargin

  def memoryRuntimeAttributeInvalid(expressionStart: Terminal) = {
    s"""ERROR: 'memory' runtime attribute should have the format "size unit" (e.g. "8 GB").
        |
        |Expression starts here (line ${expressionStart.getLine}, col ${expressionStart.getColumn}):
        |
        |${pointToSource(expressionStart)}
     """.stripMargin
  }
}
