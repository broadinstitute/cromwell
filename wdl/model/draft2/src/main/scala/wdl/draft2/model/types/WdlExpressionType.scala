package wdl.draft2.model.types

import wdl.draft2.model.WdlExpression
import wom.types.WomType
import wom.values.WomValue

case object WdlExpressionType extends WomType {
  override def toDisplayString: String = "Expression"

  override def coercion(): PartialFunction[Any, WomValue] = {
    case s: String if s.startsWith("%expr:") => WdlExpression.fromString(s.replace("%expr:", ""))
  }
}
