package cromwell.binding

import cromwell.binding.AstTools.{AstNodeName, EnhancedAstNode, EnhancedAstSeq}
import cromwell.binding.command.{ParameterCommandPart, Command}
import cromwell.parser.WdlParser._

import scala.collection.JavaConverters._
import scala.language.postfixOps

object Task {
  def apply(ast: Ast, wdlSyntaxErrorFormatter: WdlSyntaxErrorFormatter): Task = {
    val name = ast.getAttribute("name").asInstanceOf[Terminal].getSourceString

    /* Examine all inputs to the task (which currently is only command parameters e.g. ${File x})
     * And ensure that any inputs that have the same name also have the exact same definition.
     * For example having a command line of `./script ${File x} ${String x}` is conflicting.
     */
    val commandParameters = ast.findAsts(AstNodeName.CommandParameter)
    val parameterNames = commandParameters.map{_.getAttribute("name").sourceString()}.toSet
    parameterNames.foreach {name =>
      val paramsWithSameName = commandParameters.filter {_.getAttribute("name").sourceString == name}
      ensureCommandParameterAstsMatch(paramsWithSameName, ast, wdlSyntaxErrorFormatter)
    }

    val commandAsts = ast.findAsts(AstNodeName.Command)
    if (commandAsts.size != 1) throw new UnsupportedOperationException("Expecting only one Command AST")
    val command = Command(commandAsts.head, wdlSyntaxErrorFormatter)
    val outputs = ast.findAsts(AstNodeName.Output) map {TaskOutput(_, wdlSyntaxErrorFormatter)}
    new Task(name, command, outputs, buildRuntimeAttributes(ast), ast)
  }

  private def ensureCommandParameterAstsMatch(paramAsts: Seq[Ast], taskAst: Ast, wdlSyntaxErrorFormatter: WdlSyntaxErrorFormatter) = {
    paramAsts.headOption.foreach {firstParamAst =>
      val sentinal = ParameterCommandPart(firstParamAst, wdlSyntaxErrorFormatter)
      paramAsts.foreach { paramAst =>
        val parsed = ParameterCommandPart(paramAst, wdlSyntaxErrorFormatter)
        if (parsed != sentinal)
          throw new SyntaxError(wdlSyntaxErrorFormatter.parametersWithSameNameMustHaveSameDefinition(
            taskAst.getAttribute("name").asInstanceOf[Terminal],
            paramAst.getAttribute("name").asInstanceOf[Terminal],
            paramAsts.head.getAttribute("name").asInstanceOf[Terminal]
          ))
      }
    }
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
