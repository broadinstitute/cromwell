package wdl4s.cwl

import wdl4s.wdl.values.WdlValue
import wdl4s.wom.CommandPart
import wdl4s.wom.expression.IoFunctionSet


case class CwlExpressionCommandPart(expr: String) extends CommandPart {
  override def instantiate(inputsMap: Map[String, WdlValue],
                            functions: IoFunctionSet,
                            valueMapper: (WdlValue) => WdlValue): String = {
    CwlWomExpression(expr).evaluateValue(inputsMap, functions).map(_.valueString).getOrElse(???)
  }
}

