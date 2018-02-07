package wdl.draft2.types

import wdl.draft2.WdlExpression
import wom.types.WomType

case object WdlExpressionType extends WomType {
  override def toDisplayString: String = "Expression"

  override protected def coercion = {
    case s: String if s.startsWith("%expr:") => WdlExpression.fromString(s.replace("%expr:", ""))
  }
}
