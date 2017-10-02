package cromwell.backend.wdl

import cromwell.backend.BackendJobDescriptor
import cromwell.core.NoIoFunctionSet
import wdl.values.WdlValue
import wom.WomEvaluatedCallInputs
import wom.expression.IoFunctionSet

import scala.util.{Success, Try}

object Command {

  /**
    * Instantiate the command for this job descriptor.
    *
    * @param jobDescriptor jobDescriptor to instantiate the command for
    * @param callEngineFunction engine functions to use to evaluate expressions inside the command
    * @param inputsPreProcessor function to be applied to the task inputs before they are used to instantiate the command
    *                   Typically this is where localization and/or file path transformation work would be done.
    *                   The return value of the function is the inputs map that will be used to resolve variables in the command line.
    * @param valueMapper function to apply, during instantiation of the command line, after a variable is resolved
    * @return
    */
  def instantiate(jobDescriptor: BackendJobDescriptor,
                  callEngineFunction: IoFunctionSet,
                  inputsPreProcessor: WomEvaluatedCallInputs => Try[WomEvaluatedCallInputs] = (i: WomEvaluatedCallInputs) => Success(i),
                  valueMapper: WdlValue => WdlValue = identity): Try[String] = {
    inputsPreProcessor(jobDescriptor.inputDeclarations) flatMap { mappedInputs =>
      jobDescriptor.call.callable.instantiateCommand(mappedInputs, NoIoFunctionSet, valueMapper)
    }
  }
}
