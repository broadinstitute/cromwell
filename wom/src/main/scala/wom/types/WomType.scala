package wom.types

import wom.WomExpressionException
import wom.values.{WomOptionalValue, WomValue}

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
  protected def coercion: PartialFunction[Any, WomValue]
  def coercionDefined(any: Any) = coercion.isDefinedAt(any)

  /**
   * Public interface for a `Try`-wrapped conversion of an input of type `Any` to
   * a `WomValue`.
   */
  final def coerceRawValue(any: Any): Try[WomValue] = {
    any match {
      case womValue: WomValue if womValue.womType == this => Success(womValue)
      case WomOptionalValue(_, Some(v)) => coerceRawValue(v)
      case womValue: WomValue if !coercion.isDefinedAt(any) => Failure(new IllegalArgumentException(
        s"No coercion defined from wom value(s) '${WomValue.takeMaxElements(womValue, 3).toWomString}' of type" +
          s" '${womValue.womType.stableName}' to '$stableName'."))
      case _ if !coercion.isDefinedAt(any) => Failure(new IllegalArgumentException(
        s"No coercion defined from '${ScalaRunTime.stringOf(any, 3)}' of type" +
          s" '${Option(any.getClass.getCanonicalName).getOrElse(any.getClass.getName)}' to '$stableName'."))
      case _ => Try(coercion(any))
    }
  }

  final def isCoerceableFrom(otherType: WomType): Boolean = otherType match {
    case WomAnyType => true
    case WomNothingType => true
    case _ => typeSpecificIsCoerceableFrom(otherType)
  }
  protected def typeSpecificIsCoerceableFrom(otherType: WomType): Boolean = otherType == this

  /**
    * Friendly name, suitable for UIs, may change over time (if you are confident all clients will support the change)
    * @return String
    */
  def friendlyName: String = stableName

  /**
    * Stable name for call cache hashing, may expose extra information to help compute equality more granularly, must never change
    * @return String
    */
  def stableName: String

  def invalid(operation: String) = Failure(new WomExpressionException(s"Type evaluation cannot determine type from expression: $operation"))
  def add(rhs: WomType): Try[WomType] = invalid(s"$this + $rhs")
  def subtract(rhs: WomType): Try[WomType] = invalid(s"$this - $rhs")
  def multiply(rhs: WomType): Try[WomType] = invalid(s"$this * $rhs")
  def divide(rhs: WomType): Try[WomType] = invalid(s"$this / $rhs")
  def mod(rhs: WomType): Try[WomType] = invalid(s"$this % $rhs")
  def equalsType(rhs: WomType): Try[WomType] =
    if(this == rhs)
      Success(WomBooleanType)
    else
      invalid(s"$this == $rhs")
  def notEquals(rhs: WomType): Try[WomType] = equalsType(rhs) map { _ => WomBooleanType}
  def lessThan(rhs: WomType): Try[WomType] = invalid(s"$this < $rhs")
  def lessThanOrEqual(rhs: WomType): Try[WomType] = (lessThan(rhs), equalsType(rhs)) match {
    case (Success(b:WomType), _) if b == WomBooleanType => Success(WomBooleanType)
    case (_, Success(b:WomType)) if b == WomBooleanType => Success(WomBooleanType)
    case (_, _) => invalid(s"$this <= $rhs")
  }
  def greaterThan(rhs: WomType): Try[WomType] = invalid(s"$this > $rhs")
  def greaterThanOrEqual(rhs: WomType): Try[WomType] = (greaterThan(rhs), equalsType(rhs)) match {
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

  def lowestCommonSubtype(types: Iterable[WomType]): WomType = types match {
    case e if e.isEmpty => WomNothingType
    case ListOfPrimitives(primitiveType) => primitiveType
    case ListOfPairs(pairType) => pairType
    case ListOfOptionals(optionalType) => optionalType
    case ListOfMaps(mapType) => mapType
    case ListOfArrays(arrayType) => arrayType
    case ListOfObject(objectType) => objectType
    case _ => WomAnyType
  }

  private object ListOfPrimitives {
    def unapply(types: Iterable[WomType]): Option[WomType] = {
      if (types.forall(_.isInstanceOf[WomPrimitiveType])) {
        firstCommonPrimitive(types)
      } else None
    }

    val coercePriority = List(
      WomStringType, WomSingleFileType, WomUnlistedDirectoryType, WomFloatType, WomIntegerType, WomBooleanType, WomObjectType
    )

    private def firstCommonPrimitive(types: Iterable[WomType]): Option[WomType] = {
      // Types in the incoming list which everything could coerce to:
      val suppliedOptions = types.filter(t => types.forall(t.isCoerceableFrom))

      // A type not in the incoming list but which everything could coerce to nonetheless:
      lazy val unsuppliedOption: Option[WomType] = coercePriority.find(p => types.forall(p.isCoerceableFrom))

      coercePriority.find { p => suppliedOptions.toList.contains(p) } orElse { unsuppliedOption }
    }
  }


  private object ListOfPairs {
    def unapply(types: Iterable[WomType]): Option[WomPairType] = {
      if (types.forall(_.isInstanceOf[WomPairType])) {
        val pairs = types.map(_.asInstanceOf[WomPairType])
        val leftType = lowestCommonSubtype(pairs.map(_.leftType))
        val rightType = lowestCommonSubtype(pairs.map(_.rightType))
        Some(WomPairType(leftType, rightType))
      } else None
    }
  }

  private object ListOfOptionals {
    def unapply(types: Iterable[WomType]): Option[WomOptionalType] = {
      val atLeastOneOptional = types.exists {
        case _: WomOptionalType => true
        case _ => false
      }

      if (atLeastOneOptional) {
        val innerTypes = types map {
          case WomOptionalType(inner) => inner
          case nonOptional => nonOptional
        }
        Some(WomOptionalType(lowestCommonSubtype(innerTypes)))
      } else {
        None
      }
    }
  }

  private object ListOfMaps {
    def unapply(types: Iterable[WomType]): Option[WomMapType] = {
      val asMapTypes = types map {
        case m: WomMapType => Some(m)
        case _ => None
      }

      if (asMapTypes.forall(_.isDefined)) {
        val maps = asMapTypes.map(_.get)
        val keyType = lowestCommonSubtype(maps.map(_.keyType))
        val valueType = lowestCommonSubtype(maps.map(_.valueType))
        Some(WomMapType(keyType, valueType))
      } else None
    }
  }

  private object ListOfObject {
    def unapply(types: Iterable[WomType]): Option[WomObjectType.type] = {
      val asObjectTypes = types map {
        case m: WomObjectTypeLike => Some(m)
        case _ => None
      }

      if (asObjectTypes.forall(_.isDefined)) {
        Some(WomObjectType)
      } else None
    }
  }

  private object ListOfArrays {
    def unapply(types: Iterable[WomType]): Option[WomArrayType] = {
      val asArrayTypes = types map {
        case a: WomArrayType => Some(a)
        case _ => None
      }

      if (asArrayTypes.forall(_.isDefined)) {
        val arrs = asArrayTypes.map(_.get)

        /*
        WomNothingType is not coercible to any other type, yet a empty array of type WomArrayType(WomNothingType) is
        compatible with every other array type. Detect this here so we don't attempt primitive coercion on a type
        that cannot be coerced.

        Equivalent to logic found in typeSpecificIsCoerceableFrom on WomArrayType.
        */
        val memberType = lowestCommonSubtype(arrs.map(_.memberType).filterNot(_.equals(WomNothingType)))
        Some(WomArrayType(memberType))
      } else None
    }
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
