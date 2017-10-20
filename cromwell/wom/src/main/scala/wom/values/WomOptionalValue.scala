package wom.values

import wom.types.{WomOptionalType, WomType}

import scala.util.Try

case class WomOptionalValue(innerType: WomType, value: Option[WomValue]) extends WomValue {
  override val womType = WomOptionalType(innerType)
  override val toWomString = value map { _.toWomString } getOrElse "null"

  override def add(rhs: WomValue): Try[WomValue] = value match {
    case Some(lhs) => lhs.add(rhs)
    case None => emptyValueFailure("+")
  }

  override def subtract(rhs: WomValue): Try[WomValue] = value match {
    case Some(lhs) => lhs.subtract(rhs)
    case None => emptyValueFailure("-")
  }
  
  override def multiply(rhs: WomValue): Try[WomValue] = value match {
    case Some(lhs) => lhs.multiply(rhs)
    case None => emptyValueFailure("*")
  }

  override def divide(rhs: WomValue): Try[WomValue] = value match {
    case Some(lhs) => lhs.divide(rhs)
    case None => emptyValueFailure("/")
  }
  
  override def mod(rhs: WomValue): Try[WomValue] = value match {
    case Some(lhs) => lhs.mod(rhs)
    case None => emptyValueFailure("%")
  }

  override def notEquals(rhs: WomValue): Try[WomBoolean] = value match {
    case Some(lhs) => lhs.notEquals(rhs)
    case None => emptyValueFailure("!=")
  }

  override def or(rhs: WomValue): Try[WomBoolean] = value match {
    case Some(lhs) => lhs.or(rhs)
    case None => emptyValueFailure("||")
  }

  override def and(rhs: WomValue): Try[WomBoolean] = value match {
    case Some(lhs) => lhs.and(rhs)
    case None => emptyValueFailure("&&")
  }

  override def not: Try[WomValue] = value match {
    case Some(lhs) => lhs.not
    case None => emptyValueFailure("!")
  }

  override def unaryPlus: Try[WomValue] = value match {
    case Some(lhs) => lhs.unaryPlus
    case None => emptyValueFailure("+")
  }

  override def unaryMinus: Try[WomValue] = value match {
    case Some(lhs) => lhs.unaryMinus
    case None => emptyValueFailure("-")
  }

  override def equals(rhs: WomValue): Try[WomBoolean] = value match {
    case Some(lhs) => lhs.equals(rhs)
    case None => emptyValueFailure("==")
  }

  override def lessThan(rhs: WomValue): Try[WomBoolean] = value match {
    case Some(lhs) => lhs.lessThan(rhs)
    case None => emptyValueFailure("<")
  }

  override def greaterThan(rhs: WomValue): Try[WomBoolean] = value match {
    case Some(lhs) => lhs.greaterThan(rhs)
    case None => emptyValueFailure(">")
  }

  override def collectAsSeq[T <: WomValue](filterFn: PartialFunction[WomValue, T]): Seq[T] = {
    value.toList flatMap { _.collectAsSeq(filterFn) }
  }

  override final lazy val valueString: String = value match {
    case Some(v) => v.valueString
    case None => ""
  }
}

object WomOptionalValue {

  def apply(value: WomValue): WomOptionalValue = Option(value) match {
    case someValue @ Some(innerValue) => WomOptionalValue(innerValue.womType, someValue)
    case None => throw new NullPointerException(s"Cannot use apply for a null WomValue")
  }

  def none(womType: WomType) = WomOptionalValue(womType, None)
}
