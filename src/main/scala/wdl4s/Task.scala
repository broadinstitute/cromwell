package wdl4s

import java.util.regex.Pattern

import wdl4s.AstTools.{AstNodeName, EnhancedAstNode}
import wdl4s.command.{CommandPart, ParameterCommandPart, StringCommandPart}
import wdl4s.expression.WdlFunctions
import wdl4s.values.WdlValue
import wdl4s.parser.WdlParser._

import scala.annotation.tailrec
import scala.collection.JavaConverters._
import scala.language.postfixOps
import scala.util.Try

object Task {
  val Ws = Pattern.compile("[\\ \\t]+")
  def apply(ast: Ast, wdlSyntaxErrorFormatter: WdlSyntaxErrorFormatter): Task = {
    val taskNameTerminal = ast.getAttribute("name").asInstanceOf[Terminal]
    val name = taskNameTerminal.sourceString
    val declarations = ast.findAsts(AstNodeName.Declaration).map(Declaration(_, "name", wdlSyntaxErrorFormatter))
    val commandAsts = ast.findAsts(AstNodeName.Command)
    val runtimeAttributes = RuntimeAttributes(ast)
    val meta = wdlSectionToStringMap(ast, AstNodeName.Meta, wdlSyntaxErrorFormatter)
    val parameterMeta = wdlSectionToStringMap(ast, AstNodeName.ParameterMeta, wdlSyntaxErrorFormatter)

    if (commandAsts.size != 1) throw new SyntaxError(wdlSyntaxErrorFormatter.expectedExactlyOneCommandSectionPerTask(taskNameTerminal))
    val commandTemplate = commandAsts.head.getAttribute("parts").asInstanceOf[AstList].asScala.toVector map {
      case x: Terminal => new StringCommandPart(x.getSourceString)
      case x: Ast => ParameterCommandPart(x, wdlSyntaxErrorFormatter)
    }

    val outputs = ast.findAsts(AstNodeName.Output) map { TaskOutput(_, declarations, wdlSyntaxErrorFormatter) }
    Task(name, declarations, commandTemplate, runtimeAttributes, meta, parameterMeta, outputs, ast)
  }

  private def wdlSectionToStringMap(taskAst: Ast, node: String, wdlSyntaxErrorFormatter: WdlSyntaxErrorFormatter): Map[String, String] = {
    taskAst.findAsts(node) match {
      case a if a.isEmpty => Map.empty[String, String]
      case a if a.size == 1 =>
        // Yes, even 'meta {}' and 'parameter_meta {}' sections have RuntimeAttribute ASTs.
        // In hindsight, this was a poor name for the AST.
        a.head.findAsts(AstNodeName.RuntimeAttribute).map({ ast =>
          val key = ast.getAttribute("key").asInstanceOf[Terminal]
          val value = ast.getAttribute("value")
          if (!value.isInstanceOf[Terminal] || value.asInstanceOf[Terminal].getTerminalStr != "string") {
            // Keys are parsed as identifiers, but values are parsed as expressions.
            // For now, only accept expressions that are strings
            throw new SyntaxError(wdlSyntaxErrorFormatter.expressionExpectedToBeString(key))
          }
          key.sourceString -> value.sourceString
        }).toMap
      case _ => throw new SyntaxError(wdlSyntaxErrorFormatter.expectedAtMostOneSectionPerTask(node, taskAst.getAttribute("name").asInstanceOf[Terminal]))
    }
  }

  def empty: Task = new Task("taskName", Seq.empty, Seq.empty, RuntimeAttributes(Map.empty[String, Seq[String]]), Map.empty, Map.empty, Seq.empty, null)
}

/**
 * Represents a `task` declaration in a WDL file
 *
 * @param name Name of the task
 * @param declarations Any declarations (e.g. String something = "hello") defined in the task
 * @param commandTemplate Sequence of command pieces, essentially a parsed command template
 * @param runtimeAttributes 'runtime' section of a file
 * @param meta 'meta' section of a task
 * @param parameterMeta - 'parameter_meta' section of a task
 * @param outputs Set of defined outputs in the `output` section of the task
 * @param ast The syntax tree from which this was built.
 */
case class Task(name: String,
                declarations: Seq[Declaration],
                commandTemplate: Seq[CommandPart],
                runtimeAttributes: RuntimeAttributes,
                meta: Map[String, String],
                parameterMeta: Map[String, String],
                outputs: Seq[TaskOutput],
                ast: Ast) extends Executable {
  import Task._

  /**
   * Inputs to this task, as locally qualified names
   *
   * @return Seq[TaskInput] where TaskInput contains the input
   *         name & type as well as any postfix quantifiers (?, +)
   */
  val inputs: Seq[TaskInput] = for (declaration <- declarations; input <- declaration.asTaskInput) yield input
  // TODO: I think TaskInput can be replaced by Declaration

  /**
   * Given a map of task-local parameter names and WdlValues, create a command String.
   *
   * Instantiating a command line is the process of taking a command in this form:
   *
   * {{{
   *   sh script.sh ${var1} -o ${var2}
   * }}}
   *
   * This command is stored as a `Seq[CommandPart]` in the `Command` class (e.g. [sh script.sh, ${var1}, -o, ${var2}]).
   * Then, given a map of variable -> value:
   *
   * {{{
   * {
   *   "var1": "foo",
   *   "var2": "bar"
   * }
   * }}}
   *
   * It calls instantiate() on each part, and passes this map. The ParameterCommandPart are the ${var1} and ${var2}
   * pieces and they lookup var1 and var2 in that map.
   *
   * The command that's returned from Command.instantiate() is:
   *
   *
   * {{{sh script.sh foo -o bar}}}
   *
   * @param parameters Parameter values
   * @return String instantiation of the command
   */
  def instantiateCommand(parameters: CallInputs, functions: WdlFunctions[WdlValue]): Try[String] = {
    Try(normalize(commandTemplate.map(_.instantiate(declarations, parameters, functions)).mkString("")))
  }

  def commandTemplateString: String = normalize(commandTemplate.map(_.toString).mkString)

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
    parts.map(_.substring(indent)).mkString("\n")
  }

  private def leadingWhitespaceCount(s: String): Int = {
    val matcher = Ws.matcher(s)
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

  override def toString: String = s"[Task name=$name commandTemplate=$commandTemplate}]"
}
