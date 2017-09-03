package cromwell.backend.wdl

import cromwell.backend.BackendJobDescriptor
import cromwell.core.JobOutput
import wdl4s.wdl.LocallyQualifiedName
import wdl4s.wdl.values.WdlValue
import wdl4s.wom.expression.IoFunctionSet

import scala.language.postfixOps
import scala.util.{Success, Try}

object OutputEvaluator {
  def evaluateOutputs(jobDescriptor: BackendJobDescriptor,
                      ioFunctions: IoFunctionSet,
                      postMapper: WdlValue => Try[WdlValue] = v => Success(v)): Try[Map[LocallyQualifiedName, JobOutput]] = Try {
    val knownValues = jobDescriptor.inputDeclarations map {
      case (declaration, value) => declaration.fullyQualifiedName -> value
    }
    
    jobDescriptor.call.callable.outputs map { output =>
      // TODO WOM: Should evaluateValue return a Future ?
      // TODO WOM: Aggregate failures
      output.name -> JobOutput(output.expression.evaluateValue(knownValues, ioFunctions).getOrElse(throw new Exception("Output evaluation failed")))
    } toMap
  }
}
