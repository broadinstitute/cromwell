package wdl

import java.util.regex.Pattern

import cats.data.Validated.Valid
import cats.implicits._
import common.validation.ErrorOr.ErrorOr
import wdl.AstTools._
import wdl.command._
import wdl.expression.WdlFunctions
import wdl.util.StringUtil
import wdl4s.parser.WdlParser._
import wom.InstantiatedCommand
import wom.callable.Callable.{InputDefinitionWithDefault, OptionalInputDefinition, RequiredInputDefinition}
import wom.callable.{Callable, CallableTaskDefinition, TaskDefinition}
import wom.graph.LocalName
import wom.types.WomOptionalType
import wom.values.WomValue

import scala.collection.JavaConverters._
import scala.language.postfixOps

object WdlTask {
  val Ws = Pattern.compile("[\\ \\t]+")
  private implicit val instantiatedCommandMonoid = cats.derive.monoid[InstantiatedCommand]

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
    val runtimeAttributes = WdlRuntimeAttributes(ast)
    val meta = AstTools.wdlSectionToStringMap(ast, AstNodeName.Meta, wdlSyntaxErrorFormatter)
    val parameterMeta = AstTools.wdlSectionToStringMap(ast, AstNodeName.ParameterMeta, wdlSyntaxErrorFormatter)

    if (commandAsts.size != 1) throw new SyntaxError(wdlSyntaxErrorFormatter.expectedExactlyOneCommandSectionPerTask(taskNameTerminal))
    val commandTemplate = commandAsts.head.getAttribute("parts").asInstanceOf[AstList].asScala.toVector map {
      case x: Terminal => StringCommandPart(x.getSourceString)
      case x: Ast => ParameterCommandPart(x, wdlSyntaxErrorFormatter)
    }

    WdlTask(name, commandTemplate, runtimeAttributes, meta, parameterMeta, ast)
  }

  def empty: WdlTask = new WdlTask("taskName", Seq.empty, WdlRuntimeAttributes(Map.empty[String, WdlExpression]), Map.empty, Map.empty, null)

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
                   commandTemplate: Seq[WdlCommandPart],
                   runtimeAttributes: WdlRuntimeAttributes,
                   meta: Map[String, String],
                   parameterMeta: Map[String, String],
                   ast: Ast) extends WdlCallable {

  override lazy val womDefinition = Valid(buildWomTaskDefinition)

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
    * @param taskInputs  Map[String, WomValue] of inputs to this call, keys should be declarations
    * @param functions   Implementation of the WDL standard library functions to evaluate functions in expressions
    * @param valueMapper Optional WomValue => WomValue function that is called on the result of each expression
    *                    evaluation (i.e. evaluation of $${...} blocks).
    * @return String instantiation of the command
    *
    * TODO WOM this method should go away in the #2944 purge
    */
  def instantiateCommand(taskInputs: EvaluatedTaskInputs,
                         functions: WdlFunctions[WomValue],
                         valueMapper: WomValue => WomValue = identity): ErrorOr[List[InstantiatedCommand]] = {

    val mappedInputs = taskInputs.map({case (k, v) => k.unqualifiedName -> v})
    // `foldMap`: `map` over the elements of the `List[WdlCommandPart]`s, transforming each `WdlCommandPart` to an
    // `ErrorOr[InstantiatedCommand]`. Then fold the resulting `List[ErrorOr[InstantiatedCommand]]` into a single
    // `ErrorOr[InstantiatedCommand]`.
    import WdlTask.instantiatedCommandMonoid
    val fullInstantiatedCommand: ErrorOr[InstantiatedCommand] = commandTemplate.toList
      .flatTraverse(_.instantiate(declarations, mappedInputs, functions, valueMapper)).map(_.combineAll)
    
    // `normalize` the instantiation (i.e. don't break Python code indentation)
    fullInstantiatedCommand map { c => List(c.copy(commandString = StringUtil.normalize(c.commandString)))}
  }

  def commandTemplateString: String = StringUtil.normalize(commandTemplate.map(_.toString).mkString)

  override def toString: String = s"[Task name=$name commandTemplate=$commandTemplate}]"

  /**
    * Assign declaration values from the given input map.
    * Fqn must be task declaration fqns
    * e.g.:
    * task t {
    *   String s
    * }
    * inputMap = Map("t.s" -> WdlString("hello"))
    */
  // TODO WOM: Unused except in Specs
  def inputsFromMap(inputs: Map[FullyQualifiedName, WomValue]): EvaluatedTaskInputs = {
    declarations flatMap { declaration =>
      inputs collectFirst {
        case (fqn, value) if fqn == declaration.fullyQualifiedName => declaration -> value }
    } toMap
  }

  private def buildWomTaskDefinition: TaskDefinition = CallableTaskDefinition(
    name = name,
    commandTemplateBuilder = Function.const(commandTemplate.validNel),
    runtimeAttributes = runtimeAttributes.toWomRuntimeAttributes(this),
    meta = meta,
    parameterMeta = parameterMeta,
    outputs = outputs.map(_.womOutputDefinition).toList,
    inputs = buildWomInputs,
    adHocFileCreation = Set.empty,
    environmentExpressions = Map.empty
  )

  private def buildWomInputs: List[Callable.InputDefinition] = declarations collect {
    case d if d.expression.isEmpty && !d.womType.isInstanceOf[WomOptionalType] =>
      RequiredInputDefinition(LocalName(d.unqualifiedName), d.womType)
    case d if d.expression.isEmpty && d.womType.isInstanceOf[WomOptionalType] =>
      OptionalInputDefinition(LocalName(d.unqualifiedName), d.womType.asInstanceOf[WomOptionalType])
    case d if d.expression.nonEmpty =>
      InputDefinitionWithDefault(LocalName(d.unqualifiedName), d.womType, WdlWomExpression(d.expression.get, this))
  } toList
}
