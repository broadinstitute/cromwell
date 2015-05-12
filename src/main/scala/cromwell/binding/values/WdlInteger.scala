package cromwell.binding.values

import cromwell.binding.types.WdlIntegerType

case class WdlInteger(value: Integer) extends WdlPrimitive {
  val wdlType = WdlIntegerType
  def add(rhs: WdlValue): WdlValue = {
    rhs match {
      case r:WdlInteger => WdlInteger(value + r.value)
      case r:WdlString => WdlString(value + r.value)
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
