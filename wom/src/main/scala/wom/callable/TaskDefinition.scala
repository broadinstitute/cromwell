package wom.callable

import cats.data.OptionT
import cats.implicits._
import common.Checked
import common.validation.ErrorOr.ErrorOr
import common.validation.Validation._
import wdl.util.StringUtil
import wom.callable.TaskDefinition.OutputEvaluationFunction
import wom.core._
import wom.expression.{IoFunctionSet, WomExpression}
import wom.graph.GraphNodePort.OutputPort
import wom.graph.{Graph, TaskCall}
import wom.values.{WomEvaluatedCallInputs, WomGlobFile, WomValue}
import wom.{CommandPart, InstantiatedCommand, RuntimeAttributes}

import scala.concurrent.{ExecutionContext, Future}
import scala.util.Try

object TaskDefinition {

  type EvaluatedOutputs = Checked[Map[OutputPort, WomValue]]

  /*
    * Result type of the OutputEvaluationFunction. Equivalent to Future[Option[EvaluatedOutputs]]
    * - Future because it might involve I/O
    * - Option because it might not return any value (dependent on the task)
    * - Checked because the evaluation could fail (invalid types, missing values etc...)
  */
  type OutputFunctionResponse = OptionT[Future, EvaluatedOutputs]

  // Function definition to evaluate task outputs
  type OutputEvaluationFunction = (Set[OutputPort], Map[String, WomValue], IoFunctionSet, ExecutionContext) => OutputFunctionResponse
  
  object OutputEvaluationFunction {
    val none: OutputEvaluationFunction = { case (_ ,_, _, _) => OptionT[Future, EvaluatedOutputs](Future.successful(None)) }
  }

  private implicit val instantiatedCommandMonoid = cats.derive.monoid[InstantiatedCommand]
  object CommandTemplateBuilder {
    def fromValues(values: Seq[CommandPart]) = new CommandTemplateBuilder {
      override def build(inputs: WomEvaluatedCallInputs): ErrorOr[Seq[CommandPart]] = values.validNel
    }
  }
  abstract class CommandTemplateBuilder {
    def build(inputs: WomEvaluatedCallInputs): ErrorOr[Seq[CommandPart]]
  }
}

sealed trait TaskDefinition extends Callable {

  def commandTemplateBuilder: WomEvaluatedCallInputs => ErrorOr[Seq[CommandPart]]
  // TODO ErrorOrify this ? Throw for now
  def commandTemplate(taskInputs: WomEvaluatedCallInputs): Seq[CommandPart] = commandTemplateBuilder(taskInputs).toTry("Failed to build command").get
  def runtimeAttributes: RuntimeAttributes
  def meta: Map[String, String]
  def parameterMeta: Map[String, String]
  def prefixSeparator: String
  def commandPartSeparator: String
  def stdinRedirection: Option[WomExpression]
  def stdoutRedirection: Option[String]
  def stderrRedirection: Option[String]
  def adHocFileCreation: Set[WomExpression]
  def environmentExpressions: Map[String, WomExpression]
  def additionalGlob: Option[WomGlobFile]
  /**
    * Provides a custom way to evaluate outputs of the task definition.
    * Return None to leave the evaluation method to the engine.
    */
  private [wom] def customizedOutputEvaluation: OutputEvaluationFunction

  lazy val unqualifiedName: LocallyQualifiedName = name

  def instantiateCommand(taskInputs: WomEvaluatedCallInputs,
                         functions: IoFunctionSet,
                         valueMapper: WomValue => WomValue,
                         runtimeEnvironment: RuntimeEnvironment): ErrorOr[InstantiatedCommand] = {

    val mappedInputs = taskInputs.map({case (k, v) => k.localName -> v})
    import TaskDefinition.instantiatedCommandMonoid

    // Just raw command parts, no separators.
    val rawCommandParts: List[ErrorOr[InstantiatedCommand]] =
      commandTemplate(taskInputs).toList.flatMap({ commandPart =>
        commandPart.instantiate(mappedInputs, functions, valueMapper, runtimeEnvironment).sequence
      })

    // Add separator command parts and monoid smash down to one `ErrorOr[InstantiatedCommand]`.
    val instantiatedCommand: ErrorOr[InstantiatedCommand] =
      rawCommandParts.intercalate(InstantiatedCommand(commandPartSeparator).validNel)

    // `normalize` the instantiation (i.e. don't break Python code indentation)
    instantiatedCommand map { c => c.copy(commandString = StringUtil.normalize(c.commandString))}
  }

  def commandTemplateString(taskInputs: WomEvaluatedCallInputs): String = StringUtil.normalize(commandTemplate(taskInputs).map(_.toString).mkString)

  override def toString: String = {
    val template = Try(commandTemplate(Map.empty).toString()).getOrElse("Could not generate command template without inputs")
    s"[Task name=$name commandTemplate=$template]"
  }
}

/**
  * A task definition only.
  * Can be called but cannot be used in an Executable as a standalone execution.
  */
final case class CallableTaskDefinition(name: String,
                                        commandTemplateBuilder: WomEvaluatedCallInputs => ErrorOr[Seq[CommandPart]],
                                        runtimeAttributes: RuntimeAttributes,
                                        meta: Map[String, String],
                                        parameterMeta: Map[String, String],
                                        outputs: List[Callable.OutputDefinition],
                                        inputs: List[_ <: Callable.InputDefinition],
                                        adHocFileCreation: Set[WomExpression],
                                        environmentExpressions: Map[String, WomExpression],
                                        prefixSeparator: String = ".",
                                        commandPartSeparator: String = "",
                                        stdinRedirection: Option[WomExpression] = None,
                                        stdoutRedirection: Option[String] = None,
                                        stderrRedirection: Option[String] = None,
                                        additionalGlob: Option[WomGlobFile] = None,
                                        private [wom] val customizedOutputEvaluation: OutputEvaluationFunction = OutputEvaluationFunction.none
                                       ) extends TaskDefinition

/**
  * A task definition with an embedded graph.
  * Can be called from a workflow but can also be run as a standalone execution.
  */
final case class ExecutableTaskDefinition private (callableTaskDefinition: CallableTaskDefinition,
                                                   override val graph: Graph
                                                  ) extends TaskDefinition with ExecutableCallable {
  override def name = callableTaskDefinition.name
  override def inputs = callableTaskDefinition.inputs
  override def outputs = callableTaskDefinition.outputs

  override def commandTemplateBuilder = callableTaskDefinition.commandTemplateBuilder
  override def runtimeAttributes = callableTaskDefinition.runtimeAttributes
  override def meta = callableTaskDefinition.meta
  override def parameterMeta = callableTaskDefinition.parameterMeta
  override def prefixSeparator = callableTaskDefinition.prefixSeparator
  override def commandPartSeparator = callableTaskDefinition.commandPartSeparator
  override def stdinRedirection = callableTaskDefinition.stdinRedirection
  override def stdoutRedirection = callableTaskDefinition.stdoutRedirection
  override def stderrRedirection = callableTaskDefinition.stderrRedirection
  override def adHocFileCreation = callableTaskDefinition.adHocFileCreation
  override def environmentExpressions = callableTaskDefinition.environmentExpressions
  override def additionalGlob = callableTaskDefinition.additionalGlob
  override private [wom]  def customizedOutputEvaluation = callableTaskDefinition.customizedOutputEvaluation
}

object ExecutableTaskDefinition {
  def tryApply(callableTaskDefinition: CallableTaskDefinition): ErrorOr[ExecutableTaskDefinition] =
    TaskCall.graphFromDefinition(callableTaskDefinition) map { ExecutableTaskDefinition(callableTaskDefinition, _) }
}
