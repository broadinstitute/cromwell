package wom.callable

import wom.core._
import common.validation.ErrorOr.ErrorOr
import wdl.util.StringUtil
import wom.expression.IoFunctionSet
import wom.graph.{Graph, TaskCall}
import wom.values.{WomEvaluatedCallInputs, WomValue}
import wom.{CommandPart, RuntimeAttributes}

import scala.util.Try

sealed trait TaskDefinition extends Callable {

  def commandTemplate: Seq[CommandPart]
  def runtimeAttributes: RuntimeAttributes
  def meta: Map[String, String]
  def parameterMeta: Map[String, String]
  def prefixSeparator: String
  def commandPartSeparator: String
  def stdoutRedirection: Option[String]
  def stderrRedirection: Option[String]

  lazy val unqualifiedName: LocallyQualifiedName = name

  def instantiateCommand(taskInputs: WomEvaluatedCallInputs,
                         functions: IoFunctionSet,
                         valueMapper: WomValue => WomValue = identity[WomValue],
                         separate: Boolean = false): Try[String] = {
    val mappedInputs = taskInputs.map({case (k, v) => k.localName -> v})
    Try(StringUtil.normalize(commandTemplate.map(_.instantiate(mappedInputs, functions, valueMapper)).mkString(commandPartSeparator)))
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

  def tryApply(name: String,
               commandTemplate: Seq[CommandPart],
               runtimeAttributes: RuntimeAttributes,
               meta: Map[String, String],
               parameterMeta: Map[String, String],
               outputs: List[Callable.OutputDefinition],
               inputs: List[_ <: Callable.InputDefinition],
               prefixSeparator: String = ".",
               commandPartSeparator: String = "",
               stdoutRedirection: Option[String] = None,
               stderrRedirection: Option[String] = None): ErrorOr[ExecutableTaskDefinition] =
    tryApply(CallableTaskDefinition(name, commandTemplate, runtimeAttributes, meta, parameterMeta, outputs, inputs, prefixSeparator, commandPartSeparator, stdoutRedirection, stderrRedirection))
}
