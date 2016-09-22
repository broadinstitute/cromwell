package cromwell.backend

import cromwell.core.{JobOutputs, JobOutput}
import wdl4s._
import wdl4s.expression.WdlStandardLibraryFunctions
import wdl4s.util.TryUtil
import wdl4s.values.WdlValue

import scala.util.{Success, Try}

object OutputEvaluator {
  def evaluateOutputs(jobDescriptor: BackendJobDescriptor,
                      wdlFunctions: WdlStandardLibraryFunctions,
                      postMapper: WdlValue => Try[WdlValue] = v => Success(v)): Try[JobOutputs] = {
    val inputs = jobDescriptor.inputs
    val evaluatedOutputs = jobDescriptor.call.task.outputs.
      foldLeft(Map.empty[LocallyQualifiedName, Try[JobOutput]])((outputMap, output) => {
        val currentOutputs = outputMap collect {
          case (name, value) if value.isSuccess => name -> value.get.wdlValue
        }
        def lookup = (currentOutputs ++ inputs).apply _
        val coerced = output.requiredExpression.evaluate(lookup, wdlFunctions) flatMap output.wdlType.coerceRawValue
        val jobOutput = output.name -> (coerced flatMap postMapper map JobOutput)

        outputMap + jobOutput

      })

    TryUtil.sequenceMap(evaluatedOutputs, s"Workflow ${jobDescriptor.workflowDescriptor.id} post processing failed.")
  }
}
