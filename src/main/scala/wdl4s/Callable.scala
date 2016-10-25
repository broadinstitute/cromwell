package wdl4s

import wdl4s.expression.WdlStandardLibraryFunctions
import wdl4s.util.TryUtil
import wdl4s.values.WdlValue

import scala.util.{Success, Try}

trait Callable extends Scope {
  
  def unqualifiedName: String
  
  // Assumes that this will not be accessed before the children for the task are set, otherwise it will be empty
  // If that assumption proves false, make it a def or a var that is set after children are.
  lazy val outputs: Seq[TaskOutput] = children collect { case output: TaskOutput => output }

  def evaluateOutputs(inputs: EvaluatedTaskInputs,
                      wdlFunctions: WdlStandardLibraryFunctions,
                      postMapper: WdlValue => Try[WdlValue] = v => Success(v)): Try[Map[LocallyQualifiedName, WdlValue]] = {
    val fqnInputs = inputs map { case (d, v) => d.fullyQualifiedName -> v }
    val evaluatedOutputs = outputs.foldLeft(Map.empty[Scope, Try[WdlValue]])((outputMap, output) => {
      val currentOutputs = outputMap collect {
        case (outputName, value) if value.isSuccess => outputName.fullyQualifiedName -> value.get
      }
      def knownValues = currentOutputs ++ fqnInputs
      val lookup = lookupFunction(knownValues, wdlFunctions, relativeTo = output)
      val coerced = output.requiredExpression.evaluate(lookup, wdlFunctions) flatMap output.wdlType.coerceRawValue
      val jobOutput = output -> (coerced flatMap postMapper)

      outputMap + jobOutput

    }) map { case (k, v) => k.unqualifiedName -> v }

    TryUtil.sequenceMap(evaluatedOutputs, s"Failed to evaluate outputs.")
  }

}
