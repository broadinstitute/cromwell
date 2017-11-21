package cromwell.backend.wdl

import common.validation.ErrorOr._
import common.validation.Validation._
import cromwell.backend.BackendJobDescriptor
import wom.InstantiatedCommand
import wom.callable.RuntimeEnvironment
import wom.expression.IoFunctionSet
import wom.values.{WomEvaluatedCallInputs, WomValue}

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
                  valueMapper: WomValue => WomValue,
                  runtimeEnvironment: RuntimeEnvironment): ErrorOr[InstantiatedCommand] = {
    inputsPreProcessor(jobDescriptor.evaluatedTaskInputs).toErrorOr flatMap { mappedInputs =>
      jobDescriptor.taskCall.callable.instantiateCommand(mappedInputs, callEngineFunction, valueMapper, runtimeEnvironment)
    }
  }
}
