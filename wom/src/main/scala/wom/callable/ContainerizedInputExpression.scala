package wom.callable

import common.validation.ErrorOr.ErrorOr
import wom.expression.IoFunctionSet
import wom.values.WomValue

/**
  * This is an expression that uses "containerized" (a.k.a. "mapped" input values) in its evaluation.
  */
trait ContainerizedInputExpression {
  def evaluate(hostInputValues: Map[String, WomValue], containerizedInputValues: Map[String, WomValue], ioFunctionSet: IoFunctionSet): ErrorOr[WomValue]
}
