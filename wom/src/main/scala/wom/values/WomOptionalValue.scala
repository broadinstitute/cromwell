package wom.values

import cats.Applicative
import cats.syntax.functor._
import common.validation.IOChecked.IOChecked
import wom.expression.IoFunctionSet
import wom.types.{WomOptionalType, WomType}

import scala.annotation.tailrec
import scala.language.higherKinds
import scala.util.{Success, Try}

final case class WomOptionalValue(innerType: WomType, value: Option[WomValue]) extends WomValue {
  override val womType = WomOptionalType(innerType)
  override def toWomString = value map { _.toWomString } getOrElse "null"

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

  /**
    * Unpack a nested option down to a single layer of optionality
    * eg Bring Int?[...]? down to Int?
    */
  @tailrec
  def flattenOptional: WomOptionalValue = this match {
    case WomOptionalValue(_: WomOptionalType, Some(innerOptionalValue: WomOptionalValue)) => innerOptionalValue.flattenOptional
    case WomOptionalValue(innerType: WomOptionalType, None) => WomOptionalValue(innerType.baseMemberType, None)
    case _ => this
  }

  /**
    * Box up a Some optional until it reaches a certain level (eg boxing Int? up to Int???)
    */
  @tailrec
  private def boxUntilType(targetType: WomOptionalType): WomOptionalValue = {
    assert(value.isDefined)

    assert(
      targetType.baseMemberTypeIsCompatibleWith(womType.baseMemberType),
      s"base member type ${targetType.baseMemberType} and womtype ${womType.baseMemberType} are not compatible")
    if (womType.depth.equals(targetType.depth)) {
      this
    } else {
      WomOptionalValue(womType, Some(this)).boxUntilType(targetType)
    }
  }

  /**
    * Flattens this optional value down to its base form (eg Int???? goes to Int?) and then:
    * If it's a Some:
    *  - Coerce the member element to the required type (eg Int to String)
    *  - Re-box the optional to the required final format (eg String? back to String??)
    * If it's a None:
    *  - Make an unnested None with the appropriate final type (eg use String???(None) rather than String???(Some(Some(None))))
    * @param womOptionalType The final type we want to create
    */
  def coerceAndSetNestingLevel(womOptionalType: WomOptionalType) = this.flattenOptional match {
    case WomOptionalValue(_, Some(v)) =>
      womOptionalType.baseMemberType.coerceRawValue(v).map(WomOptionalValue(_).boxUntilType(womOptionalType))
    case WomOptionalValue(_, None) => Success(WomOptionalValue(womOptionalType.memberType, None))
  }

  override lazy val valueString: String = value match {
    case Some(v) => v.valueString
    case None => ""
  }
  
  def traverse[A <: WomValue, G[_]](f: WomValue => G[A])(implicit applicative: Applicative[G]) = value map { v =>
    applicative.map(f(v)) {
      WomOptionalValue(_)
    }
  } getOrElse applicative.pure(this)

  override def initialize(ioFunctionSet: IoFunctionSet): IOChecked[WomValue] = traverse(_.initialize(ioFunctionSet)).widen
}

object WomOptionalValue {

  def apply(value: WomValue): WomOptionalValue = Option(value) match {
    case someValue @ Some(innerValue) => WomOptionalValue(innerValue.womType, someValue)
    case None => throw new NullPointerException(s"Cannot use apply for a null WomValue")
  }

  def none(womType: WomType) = WomOptionalValue(womType, None)

  object Flattened {
    def unapply(value: WomValue): Option[Option[WomValue]] = value match {
      case opt: WomOptionalValue => Some(opt.flattenOptional.value)
      case _ => None
    }
  }
}
