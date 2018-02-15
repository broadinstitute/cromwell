package wdl.model.draft3.elements

import wom.values.WomPrimitive

sealed trait ExpressionElement

case class PrimitiveLiteralExpression(value: WomPrimitive) extends ExpressionElement
