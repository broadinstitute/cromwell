package wdl4s

import java.util.regex.Pattern

import wdl4s.AstTools.{AstNodeName, EnhancedAstNode}
import wdl4s.command.{CommandPart, ParameterCommandPart, StringCommandPart}
import wdl4s.expression.{NoFunctions, WdlStandardLibraryFunctionsType, WdlFunctions}
import wdl4s.parser.WdlParser._
import wdl4s.types.{WdlIntegerType, WdlType}
import wdl4s.values.WdlValue

import scala.annotation.tailrec
import scala.collection.JavaConverters._
import scala.language.postfixOps
import scala.util.{Failure, Success, Try}

object Task {
  val Ws = Pattern.compile("[\\ \\t]+")

  /** The function validateDeclaration() and the DeclarationAccumulator class are used
    * to accumulate errors and keep track of which Declarations/TaskOutputs have been examined.
    *
    * We're using this approach instead of a scalaz ValidationNel because we still want to
    * accumulate Declarations even if there was an error with that particular
    * Declaration
    */
  case class DeclarationAccumulator(errors: Seq[String] = Seq.empty, declarations: Seq[Declaration] = Seq.empty)

  def apply(ast: Ast, wdlSyntaxErrorFormatter: WdlSyntaxErrorFormatter): Task = {
    val taskNameTerminal = ast.getAttribute("name").asInstanceOf[Terminal]
    val name = taskNameTerminal.sourceString
    val declarations = ast.findAsts(AstNodeName.Declaration).map(Declaration(_, wdlSyntaxErrorFormatter))
    val commandAsts = ast.findAsts(AstNodeName.Command)
    val runtimeAttributes = RuntimeAttributes(ast)
    val meta = wdlSectionToStringMap(ast, AstNodeName.Meta, wdlSyntaxErrorFormatter)
    val parameterMeta = wdlSectionToStringMap(ast, AstNodeName.ParameterMeta, wdlSyntaxErrorFormatter)
    val outputs = ast.findAsts(AstNodeName.Output) map { TaskOutput(_, wdlSyntaxErrorFormatter) }

    if (commandAsts.size != 1) throw new SyntaxError(wdlSyntaxErrorFormatter.expectedExactlyOneCommandSectionPerTask(taskNameTerminal))
    val commandTemplate = commandAsts.head.getAttribute("parts").asInstanceOf[AstList].asScala.toVector map {
      case x: Terminal => new StringCommandPart(x.getSourceString)
      case x: Ast => ParameterCommandPart(x, wdlSyntaxErrorFormatter)
    }

    val variablesReferencedInCommand = for {
      param <- commandTemplate collect { case x: ParameterCommandPart => x }
      variable <- param.expression.variableReferences
    } yield variable

    variablesReferencedInCommand foreach { variable =>
      if (!declarations.map(_.name).contains(variable.getSourceString)) {
        throw new SyntaxError(wdlSyntaxErrorFormatter.commandExpressionContainsInvalidVariableReference(taskNameTerminal, variable))
      }
    }

    val declarationErrors = (declarations ++ outputs).foldLeft(DeclarationAccumulator())(validateDeclaration(ast, wdlSyntaxErrorFormatter))

    declarationErrors.errors match {
      case x if x.nonEmpty => throw new SyntaxError(x.mkString(s"\n${"-" * 50}\n\n"))
      case _ =>
    }

    Task(name, declarations, commandTemplate, runtimeAttributes, meta, parameterMeta, outputs, ast)
  }

  /**
    * Ensures that the current declaration doesn't have a name conflict with another declaration
    * and that the expression for the current declaration only has valid variable references in it
    *
    * @param accumulated The declarations that come lexically before 'current' as well
    *                    as the accumulated errors up until this point
    * @param current The declaration being validated
    */
  private def validateDeclaration(taskAst: Ast, wdlSyntaxErrorFormatter: WdlSyntaxErrorFormatter)
                                 (accumulated: DeclarationAccumulator, current: Declaration): DeclarationAccumulator = {
    val variableReferences = for (expr <- current.expression.toIterable; variable <- expr.variableReferences) yield variable
    val declarationAstsWithSameName = taskAst.findAsts(AstNodeName.Declaration) collect {
      case a: Ast if a.getAttribute("name").sourceString == current.name => a
    }
    val taskOutputAstsWithSameName = taskAst.findAsts(AstNodeName.Output) collect {
      case a: Ast if a.getAttribute("var").sourceString == current.name => a
    }

    val declNameTerminals = declarationAstsWithSameName.map(_.getAttribute("name").asInstanceOf[Terminal]) ++
                            taskOutputAstsWithSameName.map(_.getAttribute("var").asInstanceOf[Terminal])

    val duplicateDeclarationError = declNameTerminals match {
      case terminals if terminals.size > 1 => Option(wdlSyntaxErrorFormatter.variableDeclaredMultipleTimes(terminals(0), terminals(1)))
      case _ => None
    }

    val invalidVariableReferenceErrors = variableReferences flatMap { variable =>
      if (!accumulated.declarations.map(_.name).contains(variable.getSourceString)) {
        // .head below because we are assuming if you have a Declaration object that it must have come from a Declaration AST
        Option(wdlSyntaxErrorFormatter.declarationContainsInvalidVariableReference(
          declarationAstsWithSameName.head.getAttribute("name").asInstanceOf[Terminal],
          variable
        ))
      } else {
        None
      }
    }

    val typeErrors = (invalidVariableReferenceErrors, current.expression) match {
      case (Nil, Some(expr)) => typeCheckExpression(
        expr, current.wdlType, accumulated.declarations, declNameTerminals.head, wdlSyntaxErrorFormatter
      )
      case _ => None
    }

    DeclarationAccumulator(
      accumulated.errors ++ duplicateDeclarationError.toSeq ++ invalidVariableReferenceErrors.toSeq ++ typeErrors.toSeq,
      accumulated.declarations :+ current
    )
  }

  /** Validates that `expr`, which is assumed to come from a Declaration, is compatible with
    * `expectedType`.  If not, a string error message will be returned.
    *
    * @param expr Expression to be validated
    * @param expectedType The type ascription of the declaration
    * @param priorDeclarations Declarations that come lexically before
    * @param declNameTerminal The Terminal that represents the name of the variable in the
    *                         declaration that `expr` comes from.
    * @param wdlSyntaxErrorFormatter A syntax error formatter in case an error was found.
    * @return Some(String) if an error occurred where the String is the error message, otherwise None
    */
  private def typeCheckExpression(expr: WdlExpression, expectedType: WdlType, priorDeclarations: Seq[Declaration],
                                  declNameTerminal: Terminal, wdlSyntaxErrorFormatter: WdlSyntaxErrorFormatter): Option[String] = {

    // .get below because lookup functions should throw exceptions if they could not lookup the variable
    def lookupType(n: String) = priorDeclarations.find(_.name == n).map(_.wdlType).get

    expr.evaluateType(lookupType, new WdlStandardLibraryFunctionsType) match {
      case Success(expressionWdlType) if !expectedType.isCoerceableFrom(expressionWdlType) =>
        Option(wdlSyntaxErrorFormatter.taskOutputExpressionTypeDoesNotMatchDeclaredType(
          declNameTerminal, expressionWdlType, expectedType
        ))
      case Success(wdlType) =>
        expr.evaluate((s: String) => throw new Throwable("not implemented"), NoFunctions) match {
          case Success(value) if expectedType.coerceRawValue(value).isFailure =>
            Option(wdlSyntaxErrorFormatter.declarationExpressionNotCoerceableToTargetType(
              declNameTerminal, expectedType
            ))
          case _ => None
        }
      case Failure(ex) =>
        Option(wdlSyntaxErrorFormatter.failedToDetermineTypeOfDeclaration(declNameTerminal))
    }
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
