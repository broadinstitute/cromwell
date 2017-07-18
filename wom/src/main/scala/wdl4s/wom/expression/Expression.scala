package wdl4s.wom.expression

import wdl4s.wdl.types.WdlType
import wdl4s.wdl.values.WdlValue
import wdl4s.wom.callable.Callable.InputDefinition

sealed trait Expression {
  def inputs: Set[InputDefinition]
  def evaluate(variableValues: Map[String, WdlValue]): Unit
  def womType: WdlType
}

case class PlaceholderExpression(womType: WdlType) extends Expression {
  override def inputs: Set[InputDefinition] = ???
  override def evaluate(variableValues: Map[String, WdlValue]): Unit = ???
}
