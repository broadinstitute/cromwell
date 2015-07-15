package cromwell.binding.types

import cromwell.binding.WdlExpression

case object WdlExpressionType extends WdlType {
  override def toWdlString: String = "Expression"

  override protected def coercion = {
    case s: String if s.startsWith("%expr:") => WdlExpression.fromString(s.replace("%expr:", ""))
  }

  override def fromWdlString(rawString: String) = WdlExpression.fromString(rawString)
}
