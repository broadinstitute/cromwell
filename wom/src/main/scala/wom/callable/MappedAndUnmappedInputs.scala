package wom.callable

import common.validation.ErrorOr.ErrorOr
import wom.expression.IoFunctionSet
import wom.values.WomValue

trait MappedAndUnmappedInputs {
  def evaluateValue(inputValues: Map[String, WomValue], mappedInputValues: Map[String, WomValue], ioFunctionSet: IoFunctionSet): ErrorOr[WomValue]
}
