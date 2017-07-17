package wdl4s.wdl.values

import wdl4s.wdl.types.{WdlOptionalType, WdlType}

import scala.util.Try

case class WdlOptionalValue(innerType: WdlType, value: Option[WdlValue]) extends WdlValue {
  override val wdlType = WdlOptionalType(innerType)
  override val toWdlString = value map { _.toWdlString } getOrElse "null"

  override def add(rhs: WdlValue): Try[WdlValue] = value match {
    case Some(lhs) => lhs.add(rhs)
    case None => emptyValueFailure("+")
  }

  override def subtract(rhs: WdlValue): Try[WdlValue] = value match {
    case Some(lhs) => lhs.subtract(rhs)
    case None => emptyValueFailure("-")
  }
  
  override def multiply(rhs: WdlValue): Try[WdlValue] = value match {
    case Some(lhs) => lhs.multiply(rhs)
    case None => emptyValueFailure("*")
  }

  override def divide(rhs: WdlValue): Try[WdlValue] = value match {
    case Some(lhs) => lhs.divide(rhs)
    case None => emptyValueFailure("/")
  }
  
  override def mod(rhs: WdlValue): Try[WdlValue] = value match {
    case Some(lhs) => lhs.mod(rhs)
    case None => emptyValueFailure("%")
  }

  override def notEquals(rhs: WdlValue): Try[WdlBoolean] = value match {
    case Some(lhs) => lhs.notEquals(rhs)
    case None => emptyValueFailure("!=")
  }

  override def or(rhs: WdlValue): Try[WdlBoolean] = value match {
    case Some(lhs) => lhs.or(rhs)
    case None => emptyValueFailure("||")
  }

  override def and(rhs: WdlValue): Try[WdlBoolean] = value match {
    case Some(lhs) => lhs.and(rhs)
    case None => emptyValueFailure("&&")
  }

  override def not: Try[WdlValue] = value match {
    case Some(lhs) => lhs.not
    case None => emptyValueFailure("!")
  }

  override def unaryPlus: Try[WdlValue] = value match {
    case Some(lhs) => lhs.unaryPlus
    case None => emptyValueFailure("+")
  }

  override def unaryMinus: Try[WdlValue] = value match {
    case Some(lhs) => lhs.unaryMinus
    case None => emptyValueFailure("-")
  }

  override def equals(rhs: WdlValue): Try[WdlBoolean] = value match {
    case Some(lhs) => lhs.equals(rhs)
    case None => emptyValueFailure("==")
  }

  override def lessThan(rhs: WdlValue): Try[WdlBoolean] = value match {
    case Some(lhs) => lhs.lessThan(rhs)
    case None => emptyValueFailure("<")
  }

  override def greaterThan(rhs: WdlValue): Try[WdlBoolean] = value match {
    case Some(lhs) => lhs.greaterThan(rhs)
    case None => emptyValueFailure(">")
  }

  override def collectAsSeq[T <: WdlValue](filterFn: PartialFunction[WdlValue, T]): Seq[T] = {
    value.toList flatMap { _.collectAsSeq(filterFn) }
  }
}

object WdlOptionalValue {

  def apply(value: WdlValue): WdlOptionalValue = Option(value) match {
    case someValue @ Some(innerValue) => WdlOptionalValue(innerValue.wdlType, someValue)
    case None => throw new NullPointerException(s"Cannot use apply for a null WdlValue")
  }

  def none(wdlType: WdlType) = WdlOptionalValue(wdlType, None)
}
