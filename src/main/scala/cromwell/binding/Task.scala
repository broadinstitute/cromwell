package cromwell.binding

import cromwell.binding.AstTools.{AstNodeName, EnhancedAstNode}
import cromwell.binding.command.{Command, ParameterCommandPart}
import cromwell.parser.WdlParser._

import scala.collection.JavaConverters._
import scala.language.postfixOps

object Task {
  def apply(ast: Ast, wdlSyntaxErrorFormatter: WdlSyntaxErrorFormatter): Task = {
    val name = ast.getAttribute("name").asInstanceOf[Terminal].getSourceString
    val declarations = ast.findAsts(AstNodeName.Declaration).map(Declaration(_, "name", wdlSyntaxErrorFormatter))

    /* Examine all inputs to the task (which currently is only command parameters e.g. ${File x})
     * And ensure that any inputs that have the same name also have the exact same definition.
     * For example having a command line of `./script ${File x} ${String x}` is conflicting.
     */
    val commandParameters = ast.findAsts(AstNodeName.CommandParameter)
    val parameterNames = commandParameters.map { _.getAttribute("name").sourceString() }.toSet
    parameterNames.foreach { name =>
      val paramsWithSameName = commandParameters.filter {
        _.getAttribute("name").sourceString == name
      }
      ensureCommandParameterAstsMatch(paramsWithSameName, ast, wdlSyntaxErrorFormatter)
    }

    val commandAsts = ast.findAsts(AstNodeName.Command)
    if (commandAsts.size != 1) throw new UnsupportedOperationException("Expecting only one Command AST")
    val command = Command(commandAsts.head, wdlSyntaxErrorFormatter)
    val outputs = ast.findAsts(AstNodeName.Output) map {TaskOutput(_, wdlSyntaxErrorFormatter)}
    new Task(name, declarations, command, outputs, ast)
  }

  private def ensureCommandParameterAstsMatch(paramAsts: Seq[Ast], taskAst: Ast, wdlSyntaxErrorFormatter: WdlSyntaxErrorFormatter) = {
    paramAsts.headOption.foreach { firstParamAst =>
      val sentinel = ParameterCommandPart(firstParamAst, wdlSyntaxErrorFormatter)
      paramAsts.foreach { paramAst =>
        val parsed = ParameterCommandPart(paramAst, wdlSyntaxErrorFormatter)
        if (parsed != sentinel)
          throw new SyntaxError(wdlSyntaxErrorFormatter.parametersWithSameNameMustHaveSameDefinition(
            taskAst.getAttribute("name").asInstanceOf[Terminal],
            paramAst.getAttribute("name").asInstanceOf[Terminal],
            paramAsts.head.getAttribute("name").asInstanceOf[Terminal]
          ))
      }
    }
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
                declarations: Seq[Declaration],
                command: Command,
                outputs: Seq[TaskOutput],
                ast: Ast) extends Executable {
  /**
   * Inputs to this task, as task-local names (i.e. not fully-qualified)
   *
   * @return Map of input name to type for that input
   */
  val inputs: Seq[TaskInput] = {
    val commandInputs = command.inputs
    val declarationInputs = for(declaration <- declarations; input <- declaration.asTaskInput) yield input
    commandInputs ++ declarationInputs
  }

  val runtimeAttributes = RuntimeAttributes(ast)

  override def toString: String = s"[Task name=$name command=$command]"
}
