package wom.callable

import wom.core._
import lenthall.validation.ErrorOr.ErrorOr
import wdl.util.StringUtil
import wom.expression.IoFunctionSet
import wom.graph.{Graph, TaskCall}
import wom.values.{WomValue, WomEvaluatedCallInputs}
import wom.{CommandPart, RuntimeAttributes}

import scala.util.Try

case class TaskDefinition(name: String,
                          commandTemplate: Seq[CommandPart],
                          runtimeAttributes: RuntimeAttributes,
                          meta: Map[String, String],
                          parameterMeta: Map[String, String],
                          outputs: List[Callable.OutputDefinition],
                          inputs: List[_ <: Callable.InputDefinition],
                          prefixSeparator: String = ".",
                          commandPartSeparator: String = "") extends Callable {

  val unqualifiedName: LocallyQualifiedName = name

  override lazy val graph: ErrorOr[Graph] = TaskCall.graphFromDefinition(this)

//  def lookupFunction(knownInputs: WorkflowCoercedInputs,
//                     wdlFunctions: WdlFunctions[WomValue],
//                     outputResolver: OutputResolver = NoOutputResolver,
//                     shards: Map[Scatter, Int] = Map.empty[Scatter, Int],
//                     relativeTo: Scope = null): String => WomValue = ???

  def instantiateCommand(taskInputs: WomEvaluatedCallInputs,
                         functions: IoFunctionSet,
                         valueMapper: WomValue => WomValue = identity[WomValue],
                         separate: Boolean = false): Try[String] = {
    val mappedInputs = taskInputs.map({case (k, v) => k.localName -> v})
    // TODO: Bring back inputs: Try(StringUtil.normalize(commandTemplate.map(_.instantiate(declarations, taskInputs, functions, valueMapper)).mkString("")))
    // TODO WOM this is super WDL specific and doesn't belong here.
    Try(StringUtil.normalize(commandTemplate.map(_.instantiate(mappedInputs, functions, valueMapper)).mkString(commandPartSeparator)))
  }

  def commandTemplateString: String = StringUtil.normalize(commandTemplate.map(_.toString).mkString)

  override def toString: String = s"[Task name=$name commandTemplate=$commandTemplate}]"


  // TODO: fixup? The general version in Callable might not be good enough for Task:
//  def evaluateOutputs(inputs: EvaluatedTaskInputs,
//                      wdlFunctions: WdlStandardLibraryFunctions,
//                      postMapper: WomValue => Try[WomValue] = v => Success(v)): Try[Map[TaskOutput, WomValue]] = {
//    val fqnInputs = inputs map { case (d, v) => d.fullyQualifiedName -> v }
//    val evaluatedOutputs = outputs.foldLeft(Map.empty[TaskOutput, Try[WomValue]])((outputMap, output) => {
//      val currentOutputs = outputMap collect {
//        case (outputName, value) if value.isSuccess => outputName.fullyQualifiedName -> value.get
//      }
//      def knownValues = currentOutputs ++ fqnInputs
//      val lookup = NoLookup // TODO: Reinstate: val lookup = lookupFunction(knownValues, wdlFunctions, relativeTo = output)
//      val coerced = output.expression.evaluate(lookup, wdlFunctions) flatMap output.womType.coerceRawValue
//      val jobOutput = output -> (coerced flatMap postMapper).recoverWith {
//        case t: Throwable => Failure(new RuntimeException(s"Could not evaluate ${output.name} = ${output.expression.toWdlString}", t))
//      }
//      outputMap + jobOutput
//    }) map { case (k, v) => k -> v }
//
//    TryUtil.sequenceMap(evaluatedOutputs, "Failed to evaluate outputs.")
//  }

  /**
    * Assign declaration values from the given input map.
    * Fqn must be task declaration fqns
    * e.g.:
    * task t {
    *   String s
    * }
    * inputMap = Map("t.s" -> WdlString("hello"))
    */
//  def inputsFromMap(suppliedInputs: Map[FullyQualifiedName, WomValue]): EvaluatedTaskInputs = {
    // TODO: Reinstate:
//    inputs flatMap { i =>
//      suppliedInputs collectFirst {
//        case (fqn, value) if fqn == i.name => i -> value }
//    } toMap
//    Map.empty
//  }
}
