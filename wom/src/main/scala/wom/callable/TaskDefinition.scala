package wom.callable

import cats.implicits._
import common.validation.ErrorOr.ErrorOr
import common.validation.Validation._
import wdl.util.StringUtil
import wom.core._
import wom.expression.{IoFunctionSet, WomExpression}
import wom.graph.{Graph, TaskCall}
import wom.values.{WomEvaluatedCallInputs, WomValue}
import wom.{CommandPart, InstantiatedCommand, RuntimeAttributes}

import scala.util.Try

object TaskDefinition {
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
  def stdoutRedirection: Option[String]
  def stderrRedirection: Option[String]
  def adHocFileCreation: Set[WomExpression]
  def environmentExpressions: Map[String, WomExpression]

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
                                        stdoutRedirection: Option[String] = None,
                                        stderrRedirection: Option[String] = None
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
  override def stdoutRedirection = callableTaskDefinition.stdoutRedirection
  override def stderrRedirection = callableTaskDefinition.stderrRedirection
  override def adHocFileCreation = callableTaskDefinition.adHocFileCreation
  override def environmentExpressions = callableTaskDefinition.environmentExpressions
}

object ExecutableTaskDefinition {
  def tryApply(callableTaskDefinition: CallableTaskDefinition): ErrorOr[ExecutableTaskDefinition] =
    TaskCall.graphFromDefinition(callableTaskDefinition) map { ExecutableTaskDefinition(callableTaskDefinition, _) }
}
