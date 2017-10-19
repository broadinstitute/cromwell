package cromwell.engine.workflow.lifecycle.execution.preparation

import cats.syntax.validated._
import cromwell.core.ExecutionIndex.ExecutionIndex
import cromwell.engine.workflow.lifecycle.execution.ValueStore
import lenthall.validation.ErrorOr.ErrorOr
import shapeless.Poly1
import wom.expression.{IoFunctionSet, WomExpression}
import wom.graph.GraphNode
import wom.graph.GraphNodePort.OutputPort
import wom.values.WdlValue

object InputPointerToWdlValue extends Poly1 {
  // Function that can transform any of the coproduct types to an ErrorOr[WdlValue]
  type ToWdlValueFn = (GraphNode, Map[String, WdlValue], IoFunctionSet, ValueStore, ExecutionIndex) => ErrorOr[WdlValue]

  implicit def fromWdlValue: Case.Aux[WdlValue, ToWdlValueFn] = at[WdlValue] {
    wdlValue => (_: GraphNode, _: Map[String, WdlValue], _: IoFunctionSet, _: ValueStore, _: ExecutionIndex) =>
      wdlValue.validNel: ErrorOr[WdlValue]
  }

  implicit def fromOutputPort: Case.Aux[OutputPort, ToWdlValueFn] = at[OutputPort] {
    port => (node: GraphNode, _: Map[String, WdlValue], _: IoFunctionSet, valueStore: ValueStore, index : ExecutionIndex) =>
      // TODO WOM: This is not right, we should be able to know which one to look at
      valueStore.get(port, index)
        .orElse(valueStore.get(port, None))  
        .map(_.validNel)
        .getOrElse(s"Cannot find a value for ${port.name}".invalidNel): ErrorOr[WdlValue]
  }

  implicit def fromWomExpression: Case.Aux[WomExpression, ToWdlValueFn] = at[WomExpression] { 
    womExpression => (_: GraphNode, knownValues: Map[String, WdlValue], ioFunctions: IoFunctionSet, _: ValueStore, _ : ExecutionIndex) =>
      womExpression.evaluateValue(knownValues, ioFunctions): ErrorOr[WdlValue]
  }
}
