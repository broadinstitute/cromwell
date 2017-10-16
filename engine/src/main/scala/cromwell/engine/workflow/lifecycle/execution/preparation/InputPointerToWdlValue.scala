package cromwell.engine.workflow.lifecycle.execution.preparation

import cats.syntax.validated._
import cromwell.core.ExecutionIndex.ExecutionIndex
import cromwell.engine.workflow.lifecycle.execution.OutputStore
import lenthall.validation.ErrorOr.ErrorOr
import shapeless.Poly1
import wom.expression.{IoFunctionSet, WomExpression}
import wom.graph.GraphNodePort.OutputPort
import wom.values.WomValue

object InputPointerToWdlValue extends Poly1 {
  // Function that can transform any of the coproduct types to an ErrorOr[WomValue]
  type ToWdlValueFn = (Map[String, WomValue], IoFunctionSet, OutputStore, ExecutionIndex) => ErrorOr[WomValue]

  implicit def fromWdlValue: Case.Aux[WomValue, ToWdlValueFn] = at[WomValue] {
    womValue => (_: Map[String, WomValue], _: IoFunctionSet, _: OutputStore, _: ExecutionIndex) =>
      womValue.validNel: ErrorOr[WomValue]
  }

  implicit def fromOutputPort: Case.Aux[OutputPort, ToWdlValueFn] = at[OutputPort] {
    port => (_: Map[String, WomValue], _: IoFunctionSet, outputStore: OutputStore, index : ExecutionIndex) =>
      outputStore
        .get(port, index).map(_.validNel)
        .getOrElse(s"Cannot find a value for ${port.name}".invalidNel): ErrorOr[WomValue]
  }

  implicit def fromWomExpression: Case.Aux[WomExpression, ToWdlValueFn] = at[WomExpression] { 
    womExpression => (knownValues: Map[String, WomValue], ioFunctions: IoFunctionSet, _: OutputStore, _ : ExecutionIndex) =>
      womExpression.evaluateValue(knownValues, ioFunctions): ErrorOr[WomValue]
  }
}
