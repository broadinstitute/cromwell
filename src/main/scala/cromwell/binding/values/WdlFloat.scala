package cromwell.binding.values

import cromwell.binding.WdlExpressionException
import cromwell.binding.types.WdlFloatType

case class WdlFloat(value: Float) extends WdlPrimitive {
  val wdlType = WdlFloatType
  def add(rhs: WdlValue): WdlValue = {
    rhs match {
      case r:WdlFloat => WdlFloat(value + r.value)
      case r:WdlInteger => WdlFloat(value + r.value)
      case r:WdlString => WdlString(value + r.value)
      case _ => throw new WdlExpressionException(s"Cannot perform operation: $this + $rhs")
    }
  }
  def subtract(rhs: WdlValue): WdlValue = ???
  def multiply(rhs: WdlValue): WdlValue = ???
  def divide(rhs: WdlValue): WdlValue = ???
  def mod(rhs: WdlValue): WdlValue = ???
  def equals(rhs: WdlValue): WdlValue = ???
  def notEquals(rhs: WdlValue): WdlValue = ???
  def lessThan(rhs: WdlValue): WdlValue = ???
  def lessThanOrEqual(rhs: WdlValue): WdlValue = ???
  def greaterThan(rhs: WdlValue): WdlValue = ???
  def greaterThanOrEqual(rhs: WdlValue): WdlValue = ???
  def or(rhs: WdlValue): WdlValue = ???
  def and(rhs: WdlValue): WdlValue = ???
  def not: WdlValue = ???
  def unaryPlus: WdlValue = ???
  def unaryMinus: WdlValue = ???
  def asString = value.toString
}
