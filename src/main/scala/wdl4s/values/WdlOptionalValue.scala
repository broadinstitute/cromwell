package wdl4s.values

import wdl4s.types.{WdlOptionalType, WdlType}

import scala.util.Try

case class WdlOptionalValue(innerType: WdlType, value: Option[WdlValue]) extends WdlValue {
  override val wdlType = WdlOptionalType(innerType)
  override val toWdlString = value map { _.toWdlString } getOrElse "null"

  override def add(rhs: WdlValue): Try[WdlValue] = value match {
    case Some(lhs) => lhs.add(rhs)
    case None => emptyValue(this)
  }

  override def subtract(rhs: WdlValue): Try[WdlValue] = value match {
    case Some(lhs) => lhs.subtract(rhs)
    case None => emptyValue(this)
  }
  
  override def multiply(rhs: WdlValue): Try[WdlValue] = value match {
    case Some(lhs) => lhs.multiply(rhs)
    case None => emptyValue(this)
  }

  override def divide(rhs: WdlValue): Try[WdlValue] = value match {
    case Some(lhs) => lhs.divide(rhs)
    case None => emptyValue(this)
  }
  
  override def mod(rhs: WdlValue): Try[WdlValue] = value match {
    case Some(lhs) => lhs.mod(rhs)
    case None => emptyValue(this)
  }

  override def notEquals(rhs: WdlValue): Try[WdlBoolean] = value match {
    case Some(lhs) => lhs.notEquals(rhs)
    case None => emptyValue(this)
  }

  override def or(rhs: WdlValue): Try[WdlBoolean] = value match {
    case Some(lhs) => lhs.or(rhs)
    case None => emptyValue(this)
  }

  override def and(rhs: WdlValue): Try[WdlBoolean] = value match {
    case Some(lhs) => lhs.and(rhs)
    case None => emptyValue(this)
  }

  override def not: Try[WdlValue] = value match {
    case Some(lhs) => lhs.not
    case None => emptyValue(this)
  }

  override def unaryPlus: Try[WdlValue] = value match {
    case Some(lhs) => lhs.unaryPlus
    case None => emptyValue(this)
  }

  override def unaryMinus: Try[WdlValue] = value match {
    case Some(lhs) => lhs.unaryMinus
    case None => emptyValue(this)
  }

  override def equals(rhs: WdlValue): Try[WdlBoolean] = value match {
    case Some(lhs) => lhs.equals(rhs)
    case None => emptyValue(this)
  }

  override def lessThan(rhs: WdlValue): Try[WdlBoolean] = value match {
    case Some(lhs) => lhs.lessThan(rhs)
    case None => emptyValue(this)
  }

  override def greaterThan(rhs: WdlValue): Try[WdlBoolean] = value match {
    case Some(lhs) => lhs.greaterThan(rhs)
    case None => emptyValue(this)
  }
}

object WdlOptionalValue {

  def apply(value: WdlValue): WdlOptionalValue = Option(value) match {
    case someValue @ Some(innerValue) => WdlOptionalValue(innerValue.wdlType, someValue)
    case None => throw new NullPointerException(s"Cannot use apply for a null WdlValue")
  }

  def none(wdlType: WdlType) = WdlOptionalValue(wdlType, None)
}
