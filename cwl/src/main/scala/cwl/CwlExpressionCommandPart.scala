package cwl

import wom.CommandPart
import wom.expression.IoFunctionSet
import wom.graph.LocalName
import wom.values.WomValue


case class CwlExpressionCommandPart(expr: Expression) extends CommandPart {
  override def instantiate(inputsMap: Map[LocalName, WomValue],
                           functions: IoFunctionSet,
                           valueMapper: (WomValue) => WomValue): String = {

    val pc = ParameterContext.Empty.withInputs(inputsMap.map({ case (LocalName(localName), value) => localName -> value }), functions)

    val womValue: WomValue = expr.fold(EvaluateExpression).apply(pc)

    womValue.valueString
  }
}

