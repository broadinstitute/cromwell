package cromwell.engine.workflow.lifecycle.execution.preparation

import cats.syntax.validated._
import cromwell.core.ExecutionIndex.ExecutionIndex
import cromwell.engine.workflow.lifecycle.execution.OutputStore
import lenthall.validation.ErrorOr.ErrorOr
import shapeless.Poly1
import wom.expression.{IoFunctionSet, WomExpression}
import wom.graph.GraphNodePort.OutputPort
import wom.values.WdlValue

object InputPointerToWdlValue extends Poly1 {
  // Function that can transform any of the coproduct types to an ErrorOr[WdlValue]
  type ToWdlValueFn = (Map[String, WdlValue], IoFunctionSet, OutputStore, ExecutionIndex) => ErrorOr[WdlValue]

  implicit def fromWdlValue: Case.Aux[WdlValue, ToWdlValueFn] = at[WdlValue] {
    wdlValue => (_: Map[String, WdlValue], _: IoFunctionSet, _: OutputStore, _: ExecutionIndex) =>
      wdlValue.validNel: ErrorOr[WdlValue]
  }

  implicit def fromOutputPort: Case.Aux[OutputPort, ToWdlValueFn] = at[OutputPort] {
    port => (_: Map[String, WdlValue], _: IoFunctionSet, outputStore: OutputStore, index : ExecutionIndex) =>
      outputStore
        .get(port, index).map(_.validNel)
        .getOrElse(s"Cannot find a value for ${port.name}".invalidNel): ErrorOr[WdlValue]
  }

  implicit def fromWomExpression: Case.Aux[WomExpression, ToWdlValueFn] = at[WomExpression] { 
    womExpression => (knownValues: Map[String, WdlValue], ioFunctions: IoFunctionSet, _: OutputStore, _ : ExecutionIndex) =>
      womExpression.evaluateValue(knownValues, ioFunctions): ErrorOr[WdlValue]
  }
}
