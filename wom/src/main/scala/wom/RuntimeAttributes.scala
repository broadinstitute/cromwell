package wom

import wom.expression.WomExpression

object RuntimeAttributes {
  def empty = new RuntimeAttributes(Map.empty)
}

case class RuntimeAttributes(attributes: Map[String, WomExpression])
