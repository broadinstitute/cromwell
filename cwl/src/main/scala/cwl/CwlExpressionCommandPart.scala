package cwl

import wdl.values.WdlValue
import wom.CommandPart
import wom.expression.IoFunctionSet
import CwlWomExpression._
import wdl4s.cwl.EvaluateExpression


case class CwlExpressionCommandPart(expr: Expression) extends CommandPart {
  override def instantiate(inputsMap: Map[String, WdlValue],
                            functions: IoFunctionSet,
                            valueMapper: (WdlValue) => WdlValue): String = {

    val pc = ParameterContext.Empty.withInputs(inputsMap, functions)

    val wdlValue: WdlValue = expr.fold(EvaluateExpression).apply(pc)

    wdlValue.valueString
  }
}

