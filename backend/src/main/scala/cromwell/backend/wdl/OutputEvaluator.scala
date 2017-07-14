package cromwell.backend.wdl

import cromwell.backend.BackendJobDescriptor
import cromwell.core.JobOutput
import wdl4s.wdl.LocallyQualifiedName
import wdl4s.wdl.expression.WdlStandardLibraryFunctions
import wdl4s.wdl.values.WdlValue

import scala.util.{Success, Try}

object OutputEvaluator {
  def evaluateOutputs(jobDescriptor: BackendJobDescriptor,
                      wdlFunctions: WdlStandardLibraryFunctions,
                      postMapper: WdlValue => Try[WdlValue] = v => Success(v)): Try[Map[LocallyQualifiedName, JobOutput]] = {
    jobDescriptor.call.task.evaluateOutputs(jobDescriptor.inputDeclarations, wdlFunctions, postMapper) map { outputs =>
      outputs map {
        case (output, value) => output.unqualifiedName -> JobOutput(value)
      } 
    }
  }
}
