package wdl.types

import wdl.WdlExpression
import wom.types.WdlType

case object WdlExpressionType extends WdlType {
  override def toWdlString: String = "Expression"

  override protected def coercion = {
    case s: String if s.startsWith("%expr:") => WdlExpression.fromString(s.replace("%expr:", ""))
  }
}
