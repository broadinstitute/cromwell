package cromwell.binding

import cromwell.binding.command.{ParameterCommandPart, StringCommandPart, Command}
import cromwell.binding.types._
import cromwell.parser.AstTools.AstNodeName
import cromwell.parser.WdlParser._
import cromwell.parser.AstTools.{EnhancedAstNode, EnhancedAstSeq}
import scala.collection.JavaConverters._
import scala.language.postfixOps

object Task {
  def apply(ast: Ast, wdlSyntaxErrorFormatter: WdlSyntaxErrorFormatter): Task = {
    val name = ast.getAttribute("name").asInstanceOf[Terminal].getSourceString

    val dupeParamAsts = ast.findAsts(AstNodeName.CommandParameter).duplicatesByName
    if (dupeParamAsts.nonEmpty) {
      throw new SyntaxError(wdlSyntaxErrorFormatter.taskHasDuplicatedInputs(ast, dupeParamAsts))
    }

    val commandAsts = ast.findAsts(AstNodeName.Command)
    if (commandAsts.size != 1) throw new UnsupportedOperationException("Expecting only one Command AST")
    val command = Command(commandAsts.head)
    val outputs = ast.findAsts(AstNodeName.Output) map {TaskOutput(_)}
    new Task(name, command, outputs, buildRuntimeAttributes(ast), ast)
  }

  // TODO/FIXME: If RuntimeAttributes turned into a real type (i.e. case class) the following crap could go into its construction
  private def buildRuntimeAttributes(ast: Ast): RuntimeAttributes = {
    val asts = ast.findAsts(AstNodeName.Runtime)
    if (asts.size > 1) throw new UnsupportedOperationException("Only one runtime block may be defined per task")
    val astList = asts.headOption map {_.getAttribute("map").asInstanceOf[AstList]}
    astList map processRuntimeAttributes getOrElse Map.empty[String, String]
  }

  private def processRuntimeAttributes(astList: AstList): RuntimeAttributes = {
    astList.asScala.toVector map {a => processRuntimeAttribute(a.asInstanceOf[Ast])} toMap
  }

  private def processRuntimeAttribute(ast: Ast): RuntimeAttribute = {
    (ast.getAttribute("key").sourceString(), ast.getAttribute("value").sourceString())
  }
}

/**
 * Represents a `task` declaration in a WDL file
 *
 * @param name Name of the task
 * @param command Abstract command defined in the `command` section
 * @param outputs Set of defined outputs in the `output` section of the task
 */
case class Task(name: String,
                command: Command,
                outputs: Seq[TaskOutput],
                runtimeAttributes: RuntimeAttributes,
                ast: Ast) extends Executable { // FIXME: I dislike bundling AST here, it's only used at WdlNamespace construction time, but for now ...
  /**
   * Inputs to this task, as task-local names (i.e. not fully-qualified)
   *
   * @return Map of input name to type for that input
   */
  val inputs: Seq[TaskInput] = command.inputs

  override def toString: String = s"[Task name=$name command=$command]"
}
