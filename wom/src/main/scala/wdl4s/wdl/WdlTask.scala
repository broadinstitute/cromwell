package wdl4s.wdl

import java.util.regex.Pattern

import lenthall.util.TryUtil
import wdl4s.parser.WdlParser._
import wdl4s.wdl.AstTools._
import wdl4s.wdl.command.{CommandPart, ParameterCommandPart, StringCommandPart}
import wdl4s.wdl.expression.{WdlFunctions, WdlStandardLibraryFunctions}
import wdl4s.wdl.types.WdlOptionalType
import wdl4s.wdl.util.StringUtil
import wdl4s.wdl.values.{WdlFile, WdlValue}
import wdl4s.wom.callable.Callable.{OptionalInputDefinition, OptionalInputDefinitionWithDefault, RequiredInputDefinition}
import wdl4s.wom.callable.{Callable, TaskDefinition}
import wdl4s.wom.expression.{Expression, PlaceholderExpression}

import scala.collection.JavaConverters._
import scala.language.postfixOps
import scala.util.{Failure, Success, Try}

object WdlTask {
  val Ws = Pattern.compile("[\\ \\t]+")

  /** The function validateDeclaration() and the DeclarationAccumulator class are used
    * to accumulate errors and keep track of which Declarations/TaskOutputs have been examined.
    *
    * We're using this approach instead of a cats ValidatedNel because we still want to
    * accumulate Declarations even if there was an error with that particular
    * Declaration
    */
  case class DeclarationAccumulator(errors: Seq[String] = Seq.empty, declarations: Seq[Declaration] = Seq.empty)

  def apply(ast: Ast, wdlSyntaxErrorFormatter: WdlSyntaxErrorFormatter): WdlTask = {
    val taskNameTerminal = ast.getAttribute("name").asInstanceOf[Terminal]
    val name = taskNameTerminal.sourceString
    val commandAsts = ast.findAsts(AstNodeName.Command)
    val runtimeAttributes = RuntimeAttributes(ast)
    val meta = AstTools.wdlSectionToStringMap(ast, AstNodeName.Meta, wdlSyntaxErrorFormatter)
    val parameterMeta = AstTools.wdlSectionToStringMap(ast, AstNodeName.ParameterMeta, wdlSyntaxErrorFormatter)

    if (commandAsts.size != 1) throw new SyntaxError(wdlSyntaxErrorFormatter.expectedExactlyOneCommandSectionPerTask(taskNameTerminal))
    val commandTemplate = commandAsts.head.getAttribute("parts").asInstanceOf[AstList].asScala.toVector map {
      case x: Terminal => StringCommandPart(x.getSourceString)
      case x: Ast => ParameterCommandPart(x, wdlSyntaxErrorFormatter)
    }

    WdlTask(name, commandTemplate, runtimeAttributes, meta, parameterMeta, ast)
  }

  def empty: WdlTask = new WdlTask("taskName", Seq.empty, RuntimeAttributes(Map.empty[String, WdlExpression]), Map.empty, Map.empty, null)

}

/**
  * Represents a `task` declaration in a WDL file
  *
  * @param name              Name of the task
  * @param commandTemplate   Sequence of command pieces, essentially a parsed command template
  * @param runtimeAttributes 'runtime' section of a file
  * @param meta              'meta' section of a task
  * @param parameterMeta     - 'parameter_meta' section of a task
  * @param ast               The syntax tree from which this was built.
  */
case class WdlTask(name: String,
                   commandTemplate: Seq[CommandPart],
                   runtimeAttributes: RuntimeAttributes,
                   meta: Map[String, String],
                   parameterMeta: Map[String, String],
                   ast: Ast) extends WdlCallable {
  
  override lazy val womDefinition = buildWomTaskDefinition

  override val unqualifiedName: LocallyQualifiedName = name

  // Assumes that this will not be accessed before the children for the task are set, otherwise it will be empty
  // If that assumption proves false, make it a def or a var that is set after children are.
  override lazy val outputs: Seq[TaskOutput] = children collect { case output: TaskOutput => output }

  /**
    * Given a map of task-local parameter names and WdlValues, create a command String.
    *
    * Instantiating a command line is the process of taking a command in this form:
    *
    * {{{
    *   sh script.sh $${var1} -o $${var2}
    * }}}
    *
    * This command is stored as a `Seq[CommandPart]` in the `Command` class
    * (e.g. [sh script.sh, $${var1}, -o, $${var2}]).
    * Then, given a map of variable -> value:
    *
    * {{{
    * {
    *   "var1": "foo",
    *   "var2": "bar"
    * }
    * }}}
    *
    * It calls instantiate() on each part, and passes this map. The ParameterCommandPart are the $${var1} and $${var2}
    * pieces and they lookup var1 and var2 in that map.
    *
    * The command that's returned from Command.instantiate() is:
    *
    * {{{sh script.sh foo -o bar}}}
    *
    * @param taskInputs  Map[String, WdlValue] of inputs to this call, keys should be declarations
    * @param functions   Implementation of the WDL standard library functions to evaluate functions in expressions
    * @param valueMapper Optional WdlValue => WdlValue function that is called on the result of each expression
    *                    evaluation (i.e. evaluation of $${...} blocks).
    * @return String instantiation of the command
    */
  def instantiateCommand(taskInputs: EvaluatedTaskInputs,
                         functions: WdlFunctions[WdlValue],
                         valueMapper: WdlValue => WdlValue = (v) => v): Try[String] = {
    Try(StringUtil.normalize(commandTemplate.map(_.instantiate(declarations, taskInputs, functions, valueMapper)).mkString("")))
  }

  def commandTemplateString: String = StringUtil.normalize(commandTemplate.map(_.toString).mkString)

  override def toString: String = s"[Task name=$name commandTemplate=$commandTemplate}]"

  /**
    * A lookup function that restricts known task values to input declarations
    */
  def inputsLookupFunction(inputs: WorkflowCoercedInputs,
                           wdlFunctions: WdlFunctions[WdlValue],
                           shards: Map[Scatter, Int] = Map.empty[Scatter, Int]): String => WdlValue = {
    outputs.toList match {
      case head :: _ => super.lookupFunction(inputs, wdlFunctions, NoOutputResolver, shards, relativeTo = head)
      case _ => super.lookupFunction(inputs, wdlFunctions, NoOutputResolver, shards, this)
    }
  }

  /**
    * Finds all write_* expressions in the command line, evaluates them,
    * and return the corresponding created file for them.
    */
  def evaluateFilesFromCommand(inputs: WorkflowCoercedInputs,
                               functions: WdlFunctions[WdlValue]) = {
    val commandExpressions = commandTemplate.collect({
      case x: ParameterCommandPart => x.expression
    })

    val writeFunctionAsts = commandExpressions.map(_.ast).flatMap(x => AstTools.findAsts(x, "FunctionCall")).collect({
      case y if y.getAttribute("name").sourceString.startsWith("write_") => y
    })

    val lookup = lookupFunction(inputs, functions)
    val evaluatedExpressionMap = writeFunctionAsts map { ast =>
      val expression = WdlExpression(ast)
      val value = expression.evaluate(lookup, functions)
      expression -> value
    } toMap

    evaluatedExpressionMap collect { 
      case (k, Success(file: WdlFile)) => k -> file
    }
  }

  /**
    * Tries to determine files that will be created by the outputs of this task.
    */
  def findOutputFiles(inputs: WorkflowCoercedInputs,
                      functions: WdlFunctions[WdlValue],
                      silenceEvaluationErrors: Boolean = true): Seq[WdlFile] = {
    val outputFiles = outputs flatMap { taskOutput =>
      val lookup = lookupFunction(inputs, functions, relativeTo = taskOutput)
      taskOutput.requiredExpression.evaluateFiles(lookup, functions, taskOutput.wdlType) match {
        case Success(wdlFiles) => wdlFiles
        case Failure(ex) => if (silenceEvaluationErrors) Seq.empty[WdlFile] else throw ex
      }
    }

    outputFiles.distinct
  }

  def evaluateOutputs(inputs: EvaluatedTaskInputs,
                      wdlFunctions: WdlStandardLibraryFunctions,
                      postMapper: WdlValue => Try[WdlValue] = v => Success(v)): Try[Map[TaskOutput, WdlValue]] = {
    val fqnInputs = inputs map { case (d, v) => d.fullyQualifiedName -> v }
    val evaluatedOutputs = outputs.foldLeft(Map.empty[TaskOutput, Try[WdlValue]])((outputMap, output) => {
      val currentOutputs = outputMap collect {
        case (outputName, value) if value.isSuccess => outputName.fullyQualifiedName -> value.get
      }
      def knownValues = currentOutputs ++ fqnInputs
      val lookup = lookupFunction(knownValues, wdlFunctions, relativeTo = output)
      val coerced = output.requiredExpression.evaluate(lookup, wdlFunctions) flatMap output.wdlType.coerceRawValue
      val jobOutput = output -> (coerced flatMap postMapper).recoverWith {
        case t: Throwable => Failure(new RuntimeException(s"Could not evaluate ${output.fullyQualifiedName} = ${output.requiredExpression.toWdlString}", t))
      }
      outputMap + jobOutput
    }) map { case (k, v) => k -> v }

    TryUtil.sequenceMap(evaluatedOutputs, "Failed to evaluate outputs.")
  }

  /**
    * Assign declaration values from the given input map.
    * Fqn must be task declaration fqns
    * e.g.:
    * task t {
    *   String s
    * }
    * inputMap = Map("t.s" -> WdlString("hello"))
    */
  def inputsFromMap(inputs: Map[FullyQualifiedName, WdlValue]): EvaluatedTaskInputs = {
    declarations flatMap { declaration =>
      inputs collectFirst {
        case (fqn, value) if fqn == declaration.fullyQualifiedName => declaration -> value }
    } toMap
  }

  private def buildWomTaskDefinition: TaskDefinition = TaskDefinition(
    name,
    commandTemplate,
    runtimeAttributes,
    meta,
    parameterMeta,
    outputs.map(_.womOutputDefinition).toSet,
    buildWomInputs,
    buildWomDeclarations
  )

  private def buildWomDeclarations: List[(String, Expression)] = declarations.toList collect {
    case d if d.expression.nonEmpty =>
      d.unqualifiedName -> PlaceholderExpression(d.wdlType)
  }

  private def buildWomInputs: Set[Callable.InputDefinition] = declarations collect {
    case d if d.expression.isEmpty && !d.wdlType.isInstanceOf[WdlOptionalType] =>
      RequiredInputDefinition(d.unqualifiedName, d.wdlType)
    case d if d.expression.isEmpty && d.wdlType.isInstanceOf[WdlOptionalType] =>
      OptionalInputDefinition(d.unqualifiedName, d.wdlType.asInstanceOf[WdlOptionalType])
    case d if d.expression.nonEmpty && d.wdlType.isInstanceOf[WdlOptionalType] =>
      OptionalInputDefinitionWithDefault(d.unqualifiedName, d.wdlType.asInstanceOf[WdlOptionalType], PlaceholderExpression(d.wdlType))
  } toSet
}
