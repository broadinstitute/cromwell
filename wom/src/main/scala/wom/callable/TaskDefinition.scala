package wom.callable

import cats.implicits._
import common.validation.ErrorOr.ErrorOr
import wdl.util.StringUtil
import wom.callable.TaskDefinition._
import wom.core._
import wom.expression.IoFunctionSet
import wom.graph.{Graph, TaskCall}
import wom.values.{WomEvaluatedCallInputs, WomValue}
import wom.{CommandPart, InstantiatedCommand, RuntimeAttributes}

object TaskDefinition {
  type CommandPartRuntimeSortFunction = (WomEvaluatedCallInputs, Seq[CommandPart]) => Seq[CommandPart]
  private implicit val instantiatedCommandMonoid = cats.derive.monoid[InstantiatedCommand]
  val DefaultCommandPartRuntimeSorting: CommandPartRuntimeSortFunction = (_, original) => original
}

sealed trait TaskDefinition extends Callable {

  def commandTemplate: Seq[CommandPart]
  def runtimeAttributes: RuntimeAttributes
  def meta: Map[String, String]
  def parameterMeta: Map[String, String]
  def prefixSeparator: String
  def commandPartSeparator: String
  def stdoutRedirection: Option[String]
  def stderrRedirection: Option[String]

  /**
    * Can be overridden to re-sort the command parts just before the command is instantiated, with access to the evaluated task inputs.
    */
  def commandTemplateRuntimeSort: CommandPartRuntimeSortFunction

  lazy val unqualifiedName: LocallyQualifiedName = name

  def instantiateCommand(taskInputs: WomEvaluatedCallInputs,
                         functions: IoFunctionSet,
                         valueMapper: WomValue => WomValue,
                         runtimeEnvironment: RuntimeEnvironment): ErrorOr[InstantiatedCommand] = {

    val mappedInputs = taskInputs.map({case (k, v) => k.localName -> v})
    import TaskDefinition.instantiatedCommandMonoid

    // Just raw command parts, no separators.
    val rawCommandParts: List[ErrorOr[InstantiatedCommand]] =
      commandTemplate.toList.map(_.instantiate(mappedInputs, functions, valueMapper, runtimeEnvironment))

    // Add separator command parts and monoid smash down to one `ErrorOr[InstantiatedCommand]`.
    val instantiatedCommand: ErrorOr[InstantiatedCommand] =
      rawCommandParts.intercalate(InstantiatedCommand(commandPartSeparator).validNel)

    // `normalize` the instantiation (i.e. don't break Python code indentation)
    instantiatedCommand map { c => c.copy(commandString = StringUtil.normalize(c.commandString))}
  }

  def commandTemplateString: String = StringUtil.normalize(commandTemplate.map(_.toString).mkString)

  override def toString: String = s"[Task name=$name commandTemplate=$commandTemplate}]"
}

/**
  * A task definition only.
  * Can be called but cannot be used in an Executable as a standalone execution.
  */
final case class CallableTaskDefinition(name: String,
                                        commandTemplate: Seq[CommandPart],
                                        runtimeAttributes: RuntimeAttributes,
                                        meta: Map[String, String],
                                        parameterMeta: Map[String, String],
                                        outputs: List[Callable.OutputDefinition],
                                        inputs: List[_ <: Callable.InputDefinition],
                                        prefixSeparator: String = ".",
                                        commandPartSeparator: String = "",
                                        stdoutRedirection: Option[String] = None,
                                        stderrRedirection: Option[String] = None,
                                        commandTemplateRuntimeSort: CommandPartRuntimeSortFunction = DefaultCommandPartRuntimeSorting
                                       ) extends TaskDefinition

/**
  * A task definition with an embedded graph.
  * Can be called from a workflow but can also be run as a standalone execution.
  */
final case class ExecutableTaskDefinition private (callableTaskDefinition: CallableTaskDefinition,
                                                   override val graph: Graph,
                                                   commandTemplateRuntimeSort: CommandPartRuntimeSortFunction = DefaultCommandPartRuntimeSorting
                                                  ) extends TaskDefinition with ExecutableCallable {
  override def name = callableTaskDefinition.name
  override def inputs = callableTaskDefinition.inputs
  override def outputs = callableTaskDefinition.outputs

  override def commandTemplate = callableTaskDefinition.commandTemplate
  override def runtimeAttributes = callableTaskDefinition.runtimeAttributes
  override def meta = callableTaskDefinition.meta
  override def parameterMeta = callableTaskDefinition.parameterMeta
  override def prefixSeparator = callableTaskDefinition.prefixSeparator
  override def commandPartSeparator = callableTaskDefinition.commandPartSeparator
  override def stdoutRedirection = callableTaskDefinition.stdoutRedirection
  override def stderrRedirection = callableTaskDefinition.stderrRedirection
}

object ExecutableTaskDefinition {
  def tryApply(callableTaskDefinition: CallableTaskDefinition): ErrorOr[ExecutableTaskDefinition] =
    TaskCall.graphFromDefinition(callableTaskDefinition) map { ExecutableTaskDefinition(callableTaskDefinition, _) }
}
