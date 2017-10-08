package cwl

import wdl.values.WdlValue
import wom.CommandPart
import wom.expression.IoFunctionSet
import CwlWomExpression._
import wdl4s.cwl.EvaluateExpression
import wom.graph.LocalName


case class CwlExpressionCommandPart(expr: Expression) extends CommandPart {
  override def instantiate(inputsMap: Map[LocalName, WdlValue],
                           functions: IoFunctionSet,
                           valueMapper: (WdlValue) => WdlValue): String = {

    val pc = ParameterContext.Empty.withInputs(inputsMap.map({ case (LocalName(localName), value) => localName -> value }), functions)

    val wdlValue: WdlValue = expr.fold(EvaluateExpression).apply(pc)

    wdlValue.valueString
  }
}

