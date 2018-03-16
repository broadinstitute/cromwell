package wom.callable

import cats.data.OptionT
import cats.implicits._
import common.Checked
import common.validation.ErrorOr.ErrorOr
import common.validation.Validation._
import wdl.util.StringUtil
import wom.callable.CommandTaskDefinition.OutputEvaluationFunction
import wom.core._
import wom.expression.{IoFunctionSet, WomExpression}
import wom.graph.GraphNodePort.OutputPort
import wom.graph.{Graph, TaskCall}
import wom.values.{WomEvaluatedCallInputs, WomGlobFile, WomValue}
import wom.{CommandPart, InstantiatedCommand, RuntimeAttributes}

import scala.concurrent.{ExecutionContext, Future}
import scala.util.Try

object CommandTaskDefinition {

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

/**
  * Interface for a TaskDefinition.
  * There are 2 types of TaskDefinition:
  *   - CommandTaskDefinition
  *   - ExpressionTaskDefinition
  */
sealed trait TaskDefinition extends Callable {
  def runtimeAttributes: RuntimeAttributes
  def meta: Map[String, String]
  def parameterMeta: Map[String, String]

  /**
    * Transform the Callable TaskDefinition to an ExecutableCallable that can be executed on its own.
    */
  def toExecutable: ErrorOr[ExecutableCallable]
  /**
    * Provides a custom way to evaluate outputs of the task definition.
    * Return None to leave the evaluation method to the engine.
    */
  private [wom] def customizedOutputEvaluation: OutputEvaluationFunction
}

/**
  * A task definition for a command line.
  * Can be Callable only or CallableExecutable
  */
sealed trait CommandTaskDefinition extends TaskDefinition {
  def stdoutOverride: Option[WomExpression]
  def stderrOverride: Option[WomExpression]
  def commandTemplateBuilder: WomEvaluatedCallInputs => ErrorOr[Seq[CommandPart]]
  // TODO ErrorOrify this ? Throw for now
  def commandTemplate(taskInputs: WomEvaluatedCallInputs): Seq[CommandPart] = commandTemplateBuilder(taskInputs).toTry("Failed to build command").get
  def prefixSeparator: String
  def commandPartSeparator: String
  def stdinRedirection: Option[WomExpression]
  def adHocFileCreation: Set[ContainerizedInputExpression]
  def environmentExpressions: Map[String, WomExpression]
  def additionalGlob: Option[WomGlobFile]
  def homeOverride: Option[RuntimeEnvironment => String]
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

    val inputsByLocalName = taskInputs map { case (k, v) => k.localName -> v }
    val valueMappedInputsByLocalName = inputsByLocalName map { case (k, v) => k -> valueMapper(v) }
    import CommandTaskDefinition.instantiatedCommandMonoid

    // Just raw command parts, no separators.
    val rawCommandParts: List[ErrorOr[InstantiatedCommand]] =
      commandTemplate(taskInputs).toList.flatMap({ commandPart =>
        commandPart.instantiate(inputsByLocalName, functions, valueMapper, runtimeEnvironment).sequence
      })

    // Add separator command parts and monoid smash down to one `ErrorOr[InstantiatedCommand]`.
    val instantiatedCommand: ErrorOr[InstantiatedCommand] =
      rawCommandParts.intercalate(InstantiatedCommand(commandPartSeparator).validNel)

    // `normalize` the instantiation (i.e. don't break Python code indentation) and add in the inputs.
    instantiatedCommand map { c => c.copy(
      commandString = StringUtil.normalize(c.commandString),
      preprocessedInputs = inputsByLocalName.toList,
      valueMappedPreprocessedInputs = valueMappedInputsByLocalName.toList
    )}
  }

  def commandTemplateString(taskInputs: WomEvaluatedCallInputs): String = StringUtil.normalize(commandTemplate(taskInputs).map(_.toString).mkString)

  override def toString: String = {
    val template = Try(commandTemplate(Map.empty).toString()).getOrElse("Could not generate command template without inputs")
    s"[Task name=$name commandTemplate=$template]"
  }
}

/**
  * A command task definition only.
  * Can be called but cannot be used in an Executable as a standalone execution.
  */
final case class CallableTaskDefinition(name: String,
                                        commandTemplateBuilder: WomEvaluatedCallInputs => ErrorOr[Seq[CommandPart]],
                                        runtimeAttributes: RuntimeAttributes,
                                        meta: Map[String, String],
                                        parameterMeta: Map[String, String],
                                        outputs: List[Callable.OutputDefinition],
                                        inputs: List[_ <: Callable.InputDefinition],
                                        adHocFileCreation: Set[ContainerizedInputExpression],
                                        environmentExpressions: Map[String, WomExpression],
                                        prefixSeparator: String = ".",
                                        commandPartSeparator: String = "",
                                        stdinRedirection: Option[WomExpression] = None,
                                        stdoutOverride: Option[WomExpression] = None,
                                        stderrOverride: Option[WomExpression] = None,
                                        additionalGlob: Option[WomGlobFile] = None,
                                        private [wom] val customizedOutputEvaluation: OutputEvaluationFunction = OutputEvaluationFunction.none,
                                        homeOverride: Option[RuntimeEnvironment => String] = None
                                       ) extends CommandTaskDefinition {
  def toExecutable: ErrorOr[ExecutableTaskDefinition] = TaskCall.graphFromDefinition(this) map { ExecutableTaskDefinition(this, _) }
}

/**
  * A command task definition with an embedded graph.
  * Can be called from a workflow but can also be run as a standalone execution.
  */
final case class ExecutableTaskDefinition private (callableTaskDefinition: CallableTaskDefinition,
                                                   override val graph: Graph
                                                  ) extends CommandTaskDefinition with ExecutableCallable {
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
  override def stdoutOverride = callableTaskDefinition.stdoutOverride
  override def stderrOverride = callableTaskDefinition.stderrOverride
  override def adHocFileCreation = callableTaskDefinition.adHocFileCreation
  override def environmentExpressions = callableTaskDefinition.environmentExpressions
  override def additionalGlob = callableTaskDefinition.additionalGlob
  override private [wom]  def customizedOutputEvaluation = callableTaskDefinition.customizedOutputEvaluation
  override def toExecutable = this.validNel
  override def homeOverride = callableTaskDefinition.homeOverride
}

sealed trait ExpressionTaskDefinition extends TaskDefinition {
  def evaluate: (Map[String, WomValue], IoFunctionSet, List[OutputPort]) => Checked[Map[OutputPort, WomValue]]
}


/**
  * An expression task definition only.
  * Can be called but cannot be used in an Executable as a standalone execution.
  */
final case class CallableExpressionTaskDefinition(name: String,
                                                  evaluate: (Map[String, WomValue], IoFunctionSet, List[OutputPort]) => Checked[Map[OutputPort, WomValue]],
                                                  runtimeAttributes: RuntimeAttributes,
                                                  meta: Map[String, String],
                                                  parameterMeta: Map[String, String],
                                                  outputs: List[Callable.OutputDefinition],
                                                  inputs: List[_ <: Callable.InputDefinition],
                                                  prefixSeparator: String = ".",
                                                  private [wom] val customizedOutputEvaluation: OutputEvaluationFunction = OutputEvaluationFunction.none
                                       ) extends ExpressionTaskDefinition {
  def toExecutable: ErrorOr[ExecutableExpressionTaskDefinition] = TaskCall.graphFromDefinition(this) map { ExecutableExpressionTaskDefinition(this, _) }
}

/**
  * An expression task definition with an embedded graph.
  * Can be called from a workflow but can also be run as a standalone execution.
  */
final case class ExecutableExpressionTaskDefinition private (callableTaskDefinition: CallableExpressionTaskDefinition,
                                                   override val graph: Graph
                                                  ) extends ExpressionTaskDefinition with ExecutableCallable {
  override def name = callableTaskDefinition.name
  override def inputs = callableTaskDefinition.inputs
  override def outputs = callableTaskDefinition.outputs

  override def evaluate = callableTaskDefinition.evaluate
  override def runtimeAttributes = callableTaskDefinition.runtimeAttributes
  override def meta = callableTaskDefinition.meta
  override def parameterMeta = callableTaskDefinition.parameterMeta
  override def toExecutable = this.validNel
  override private [wom]  def customizedOutputEvaluation = callableTaskDefinition.customizedOutputEvaluation
}
