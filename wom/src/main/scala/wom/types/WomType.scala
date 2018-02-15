package wom.types

import wom.WomExpressionException
import wom.values.WomValue

import scala.runtime.ScalaRunTime
import scala.util.{Failure, Success, Try}

class WomTypeException(message: String) extends RuntimeException(message)

trait WomType {

  /**
   * Method to be overridden by implementation classes defining a partial function
   * for the conversion of raw input values to specific implementation class value types.
   * i.e.  `WomBooleanType` should define a partial function that knows how to
   * construct `WomBoolean`s for inputs of supported types and contents.  Values for which
   * the partial function is not defined are assumed to not be convertible to the target type.
   */
  def coercion: PartialFunction[Any, WomValue]
  def coercionDefined(any: Any) = coercion.isDefinedAt(any)

  /**
   * Public interface for a `Try`-wrapped conversion of an input of type `Any` to
   * a `WomValue`.
   */
  final def coerceRawValue(any: Any): Try[WomValue] = {
    any match {
      case womValue: WomValue if womValue.womType == this => Success(womValue)
      case womValue: WomValue if !coercion.isDefinedAt(any) => Failure(new IllegalArgumentException(
        s"No coercion defined from '$womValue' of type" +
          s" '${womValue.womType.toDisplayString}' to '$toDisplayString'."))
      case _ if !coercion.isDefinedAt(any) => Failure(new IllegalArgumentException(
        s"No coercion defined from '${ScalaRunTime.stringOf(any, 3)}' of type" +
          s" '${Option(any.getClass.getCanonicalName).getOrElse(any.getClass.getName)}' to '$toDisplayString'."))
      case _ => Try(coercion(any))
    }
  }

  final def isCoerceableFrom(otherType: WomType): Boolean = otherType match {
    case WomAnyType => true
    case _ => typeSpecificIsCoerceableFrom(otherType)
  }
  protected def typeSpecificIsCoerceableFrom(otherType: WomType): Boolean = otherType == this

  def toDisplayString: String

  def invalid(operation: String) = Failure(new WomExpressionException(s"Type evaluation cannot determine type from expression: $operation"))
  def add(rhs: WomType): Try[WomType] = invalid(s"$this + $rhs")
  def subtract(rhs: WomType): Try[WomType] = invalid(s"$this - $rhs")
  def multiply(rhs: WomType): Try[WomType] = invalid(s"$this * $rhs")
  def divide(rhs: WomType): Try[WomType] = invalid(s"$this / $rhs")
  def mod(rhs: WomType): Try[WomType] = invalid(s"$this % $rhs")
  def equals(rhs: WomType): Try[WomType] = invalid(s"$this == $rhs")
  def notEquals(rhs: WomType): Try[WomType] = equals(rhs) map { _ => WomBooleanType}
  def lessThan(rhs: WomType): Try[WomType] = invalid(s"$this < $rhs")
  def lessThanOrEqual(rhs: WomType): Try[WomType] = (lessThan(rhs), equals(rhs)) match {
    case (Success(b:WomType), _) if b == WomBooleanType => Success(WomBooleanType)
    case (_, Success(b:WomType)) if b == WomBooleanType => Success(WomBooleanType)
    case (_, _) => invalid(s"$this <= $rhs")
  }
  def greaterThan(rhs: WomType): Try[WomType] = invalid(s"$this > $rhs")
  def greaterThanOrEqual(rhs: WomType): Try[WomType] = (greaterThan(rhs), equals(rhs)) match {
    case (Success(b:WomType), _) if b == WomBooleanType => Success(WomBooleanType)
    case (_, Success(b:WomType)) if b == WomBooleanType => Success(WomBooleanType)
    case (_, _) => invalid(s"$this >= $rhs")
  }
  def or(rhs: WomType): Try[WomType] = invalid(s"$this || $rhs")
  def and(rhs: WomType): Try[WomType] = invalid(s"$this && $rhs")
  def not: Try[WomType] = invalid(s"!$this")
  def unaryPlus: Try[WomType] = invalid(s"+$this")
  def unaryMinus: Try[WomType] = invalid(s"-$this")
}

object WomType {
  /* This is in the order of coercion from non-wom types */
  val womTypeCoercionOrder: Seq[WomType] = Seq(
    WomStringType, WomIntegerType, WomFloatType, WomPairType(WomAnyType, WomAnyType), WomMapType(WomAnyType, WomAnyType),
    WomArrayType(WomAnyType), WomBooleanType, WomObjectType,
    // Putting optional type last means we'll only coerce to it for JsNull.
    // That should be OK because every other type X can coerce into X? later if it needs to.
    WomOptionalType(WomAnyType)
  )

  def homogeneousTypeFromValues(values: Iterable[WomValue]): WomType =
    homogeneousTypeFromTypes(values.map(_.womType))

  def homogeneousTypeFromTypes(types: Iterable[WomType]): WomType = {
    types.toSet match {
      case s if s.isEmpty => WomNothingType
      case s if s.size == 1 => s.head
      case _ => lowestCommonSubtype(types)
    }
  }

  def lowestCommonSubtype(types: Iterable[WomType]): WomType = {
    types.collectFirst {
      case t1 if types.forall(t2 => t1.isCoerceableFrom(t2)) => t1
    } getOrElse WomAnyType
  }

  object RecursiveType {
    def unapply(in: WomType): Option[WomType] = in match {
      case WomOptionalType(other) => Some(other)
      case WomMaybeEmptyArrayType(other) => Some(other)
      case WomNonEmptyArrayType(other) => Some(other)
      case _ => None
    }
  }
}
