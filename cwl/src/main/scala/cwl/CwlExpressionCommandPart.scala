package cwl

import wdl.values.WdlValue
import wom.CommandPart
import wom.expression.IoFunctionSet


case class CwlExpressionCommandPart(expr: String) extends CommandPart {
  override def instantiate(inputsMap: Map[String, WdlValue],
                            functions: IoFunctionSet,
                            valueMapper: (WdlValue) => WdlValue): String = {
    CwlWomExpression(expr).evaluateValue(inputsMap, functions).map(_.valueString).getOrElse(???)
  }
}

