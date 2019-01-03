package wom.callable

import common.validation.ErrorOr.ErrorOr
import wom.expression.IoFunctionSet
import wom.values.{WomFile, WomValue}

/**
  * This is an expression that uses "containerized" (a.k.a. "mapped" input values) in its evaluation.
  */
trait ContainerizedInputExpression {
  def evaluate(hostInputValues: Map[String, WomValue],
               containerizedInputValues: Map[String, WomValue],
               ioFunctionSet: IoFunctionSet): ErrorOr[List[AdHocValue]]
}

final case class AdHocValue(womValue: WomFile, alternativeName: Option[String], inputName: Option[String])
