package wdl.draft2.model.types

import wdl.draft2.model.WdlExpression
import wom.types.WomType

case object WdlExpressionType extends WomType {
  override def stableName: String = "Expression"

  override protected def coercion = {
    case s: String if s.startsWith("%expr:") => WdlExpression.fromString(s.replace("%expr:", ""))
  }
}
