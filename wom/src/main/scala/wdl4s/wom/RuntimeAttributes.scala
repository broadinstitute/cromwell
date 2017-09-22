package wdl4s.wom

import wdl4s.wom.expression.WomExpression

object RuntimeAttributes {
  def empty = new RuntimeAttributes(Map.empty)
}

case class RuntimeAttributes(attributes: Map[String, WomExpression])
