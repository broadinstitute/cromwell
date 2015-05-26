package cromwell.binding

import java.io.File

import cromwell.binding.command._
import cromwell.binding.types.{WdlFileType, WdlIntegerType, WdlStringType, WdlType}
import cromwell.binding.values.WdlValue
import cromwell.parser.WdlParser
import cromwell.parser.WdlParser._
import cromwell.util.FileUtil

import scala.collection.JavaConverters._
import scala.util.{Failure, Try}

/**
 * Main interface into the `cromwell.binding` package.
 *
 * Example usage:
 *
 * {{{
 * val binding = WdlBinding.process(new File("/path/to/file.wdl"))
 * binding.workflow.calls foreach { call =>
 *      println(call)
 * }
 * }}}
 */
object WdlBinding {
  /**
   * Given a pointer to a WDL file, parse the text and build Workflow and Task
   * objects.
   *
   * @param wdlFile The file to parse/process
   * @return WdlBinding object with the parsed results
   * @throws WdlParser.SyntaxError if there was a problem parsing the source code
   * @throws UnsupportedOperationException if an error occurred constructing the
   *                                       Workflow and Task objects
   *
   */
  def process(wdlFile: File): WdlBinding = new WdlBinding(WdlBinding.getAst(wdlFile))

  def process(wdlSource: WdlSource): WdlBinding = new WdlBinding(WdlBinding.getAst(wdlSource, "string"))

  def process(wdlSource: WdlSource, resource: String): WdlBinding = new WdlBinding(WdlBinding.getAst(wdlSource, resource))

  /**
   * Given a WDL file, this will simply parse it and return the syntax tree
   * @param wdlFile The file to parse
   * @return an Abstract Syntax Tree (WdlParser.Ast) representing the structure of the code
   * @throws WdlParser.SyntaxError if there was a problem parsing the source code
   */
  def getAst(wdlFile: File): Ast = getAst(FileUtil.slurp(wdlFile), wdlFile.getName)

  def getAst(wdlSource: WdlSource, resource: String): Ast = {
    val parser = new WdlParser()
    val tokens = parser.lex(wdlSource, resource)
    val syntaxErrorFormatter = new WdlSyntaxErrorFormatter(wdlSource)
    validateAst(
      parser.parse(tokens, syntaxErrorFormatter).toAst.asInstanceOf[Ast],
      syntaxErrorFormatter
    )
  }

  def sourceString(astNode: AstNode): String = astNode.asInstanceOf[Terminal].getSourceString

  def findAsts(ast: AstNode, name: String): Set[Ast] = {
    ast match {
      case x: Ast =>
        val thisAst = if (x.getName.equals(name)) Set(x) else Set.empty[Ast]
        x.getAttributes.values().asScala.flatMap(findAsts(_, name)).toSet[Ast] ++ thisAst
      case x: AstList => x.asScala.toVector.flatMap(findAsts(_, name)).toSet[Ast]
      case x: Terminal => Set.empty[Ast]
      case _ => Set.empty[Ast]
    }
  }

  def getCallInput(ast: Ast): Map[String, WdlExpression] = getCallInputAsts(ast).map(processCallInput).toMap

  def getCallInputAsts(ast: Ast): Set[Ast] = {
    findAsts(ast, AstNodeName.Inputs) match {
      case x: Set[Ast] if x.size == 1 => WdlBinding.findAsts(x.head.getAttribute("map"), AstNodeName.IOMapping)
      case _ => Set.empty[Ast]
    }
  }

  /**
   * Validate the following things about the AST:
   *
   * 1) All `Call` blocks reference tasks that exist
   * 2) All `Call` inputs reference actual variables on the corresponding task
   * 3) Tasks do not have duplicate inputs
   * 4) `Call` input expressions (right-hand side) should only use the MemberAccess
   *    syntax (e.g: x.y) on WdlObjects (which include other `Call` invocations)
   * 5) `Call` input expressions (right-hand side) should only reference identifiers
   *    that will resolve when evaluated
   *
   * @param ast AST to validate
   */
  private def validateAst(ast: Ast, wdlSyntaxErrorFormatter: WdlSyntaxErrorFormatter): Ast = {
    val callAsts = WdlBinding.findAsts(ast, AstNodeName.Call)
    val taskAsts = WdlBinding.findAsts(ast, AstNodeName.Task)
    val workflowAsts = WdlBinding.findAsts(ast, AstNodeName.Workflow)
    if (workflowAsts.size > 1) throw new SyntaxError(wdlSyntaxErrorFormatter.tooManyWorkflows(workflowAsts.asJava))
    callAsts foreach { callAst =>
      val taskName = sourceString(callAst.getAttribute("task"))
      val taskAst = taskAsts.find { taskAst => sourceString(taskAst.getAttribute("name")) == taskName } getOrElse {
        throw new SyntaxError(wdlSyntaxErrorFormatter.callReferencesBadTaskName(callAst, taskName))
      }

      /* TODO: This is only counting inputs that are defined on the command line */
      val taskInputs = WdlBinding.findAsts(taskAst, AstNodeName.CommandParameter).map { paramAst =>
        sourceString(paramAst.getAttribute("name"))
      }
      WdlBinding.getCallInputAsts(callAst) foreach { callInput =>
        val callInputKey = sourceString(callInput.getAttribute("key"))
        if (!taskInputs.contains(callInputKey)) throw new SyntaxError(wdlSyntaxErrorFormatter.callReferencesBadTaskInput(callInput, taskAst));
      }
    }
    ast
  }

  private def getTask(ast: Ast): Task = {
    val name = ast.getAttribute("name").asInstanceOf[Terminal].getSourceString
    val commandAsts = findAsts(ast, AstNodeName.Command)
    if (commandAsts.size != 1) throw new UnsupportedOperationException("Expecting only one Command AST")
    val command = getCommand(commandAsts.head.getAttribute("parts").asInstanceOf[AstList])
    val outputs = WdlBinding.findAsts(ast, AstNodeName.Output).map(WdlBinding.getTaskOutput)
    new Task(name, command, outputs)
  }

  private def getTaskOutput(ast: Ast): TaskOutput = {
    val wdlType = getWdlType(ast.getAttribute("type"))
    val name = ast.getAttribute("var").asInstanceOf[Terminal].getSourceString
    val expression = ast.getAttribute("expression")
    new TaskOutput(name, wdlType, new WdlExpression(expression))
  }

  private def getCommand(astList: AstList): Command = {
    val parts = astList.asScala.toVector.map {
      case x: Terminal => new StringCommandPart(x.getSourceString)
      case x: Ast => getCommandParameter(x)
    }
    new Command(parts)
  }

  private def getCommandParameter(ast: Ast): ParameterCommandPart = {
    val wdlType = getWdlType(ast.getAttribute("type"))
    val name = ast.getAttribute("name").asInstanceOf[Terminal].getSourceString
    new ParameterCommandPart(wdlType, name)
  }

  private def getWdlType(ast: AstNode): WdlType = {
    ast match {
      case t: Terminal if t.getSourceString == WdlFileType.toWdlString => WdlFileType
      case t: Terminal if t.getSourceString == WdlStringType.toWdlString => WdlStringType
      case t: Terminal if t.getSourceString == WdlIntegerType.toWdlString => WdlIntegerType
      case null => WdlStringType
      case _ => throw new UnsupportedOperationException("Implement this later for compound types")
    }
  }

  private def processCallInput(ast: Ast): (String, WdlExpression) = {
    val key = ast.getAttribute("key").asInstanceOf[Terminal].getSourceString
    val expression = new WdlExpression(ast.getAttribute("value"))
    (key, expression)
  }

  object AstNodeName {
    val Task = "Task"
    val Workflow = "Workflow"
    val Command = "RawCommand"
    val Output = "Output"
    val CommandParameter = "CommandParameter"
    val Call = "Call"
    val IOMapping = "IOMapping"
    val Inputs = "Inputs"
  }
}

case class WdlBinding(ast: Ast) {

  import WdlBinding.AstNodeName

  /**
   * All `task` definitions in the WDL file
   */
  val tasks = WdlBinding.findAsts(ast, AstNodeName.Task).map(WdlBinding.getTask)

  /**
   * All `workflow` definitions in the WDL file.  There should be only one
   * which is checked in validateAst()
   */
  val workflowAsts = WdlBinding.findAsts(ast, AstNodeName.Workflow)

  /**
   * The `workflow` definition in the WDL file.  Currently, only one workflow is defined
   * per WDL file or else it's considered an error and an exception will be thrown
   */
  val workflow = processWorkflow(workflowAsts.head)

  workflow.calls.foreach {
    _.setParent(workflow)
  }

  private def processWorkflow(ast: Ast): Workflow = {
    val name = ast.getAttribute("name").asInstanceOf[Terminal].getSourceString
    val calls = WdlBinding.findAsts(ast, AstNodeName.Call).map(processCall)
    new Workflow(name, calls)
  }

  private def processCall(ast: Ast): Call = {
    val alias: Option[String] = ast.getAttribute("alias") match {
      case x: Terminal => Option(x.getSourceString)
      case _ => None
    }
    val taskName = ast.getAttribute("task").asInstanceOf[Terminal].getSourceString
    val task = findTask(taskName) getOrElse {
      throw new UnsupportedOperationException("Cannot find task with name: " + taskName)
    }
    val inputs = WdlBinding.getCallInput(ast)
    new Call(alias, task, inputs.toMap)
  }

  def findTask(name: String): Option[Task] = tasks.find(_.name == name)

  /**
   * Confirm all required inputs are present and attempt to coerce raw inputs to `WdlValue`s.
   * This can fail if required raw inputs are missing or if the values for a specified raw input
   * cannot be coerced to the target type of the input as specified in the binding.
   */
  def confirmAndCoerceRawInputs(rawInputs: Map[FullyQualifiedName, Any]): Try[Map[FullyQualifiedName, WdlValue]] = {
    val tryCoercedValues = workflow.inputs.map { case (fqn, wdlType) =>
      val tryValue = if (!rawInputs.contains(fqn)) {
        Failure(new UnsatisfiedInputsException(s"Required workflow input '$fqn' not specified."))
      } else {
        wdlType.coerceRawValue(rawInputs.get(fqn).get)
      }
      fqn -> tryValue
    }

    val (successes, failures) = tryCoercedValues.partition { case (_, tryValue) => tryValue.isSuccess }
    if (failures.isEmpty) {
      Try(successes.map { case (key, tryValue) => key -> tryValue.get })
    } else {
      val message = failures.values.collect { case f: Failure[_] => f.exception.getMessage }.mkString("\n")
      Failure(new UnsatisfiedInputsException(message))
    }
  }
}
