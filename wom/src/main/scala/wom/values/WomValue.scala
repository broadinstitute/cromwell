package wom.values

import common.collections.EnhancedCollections._
import common.validation.IOChecked
import common.validation.IOChecked.IOChecked
import wom.expression.{IoFunctionSet, ValueAsAnExpression}
import wom.types.WomType
import wom.{OptionalNotSuppliedException, WomExpressionException}

import scala.collection.immutable.TreeMap
import scala.util.{Failure, Try}

trait WomValue {
  type ConcreteType <: WomValue
  def womType: WomType
  def invalid(operation: String) = Failure(new WomExpressionException(s"Cannot perform operation: $operation"))
  def emptyValueFailure(operationName: String) = Failure(OptionalNotSuppliedException(operationName))
  def evaluateIfDefined[A <: WomValue](operationName: String, optionalValue: WomOptionalValue, operation: WomValue => Try[A]): Try[A] = optionalValue match {
    case WomOptionalValue(_, Some(v)) => operation(v)
    case _ => emptyValueFailure(operationName)
  }

  def add(rhs: WomValue): Try[WomValue] = invalid(s"$this + $rhs")
  def subtract(rhs: WomValue): Try[WomValue] = invalid(s"$this - $rhs")
  def multiply(rhs: WomValue): Try[WomValue] = invalid(s"$this * $rhs")
  def divide(rhs: WomValue): Try[WomValue] = invalid(s"$this / $rhs")
  def mod(rhs: WomValue): Try[WomValue] = invalid(s"$this % $rhs")
  def equals(rhs: WomValue): Try[WomBoolean] = invalid(s"$this == $rhs")
  def notEquals(rhs: WomValue): Try[WomBoolean] = equals(rhs).map{ x => WomBoolean(!x.value)}
  def lessThan(rhs: WomValue): Try[WomBoolean] = invalid(s"$this < $rhs")
  def lessThanOrEqual(rhs: WomValue): Try[WomBoolean] =
    Try(WomBoolean(Seq(lessThan _, equals _).exists{ p => p(rhs).get == WomBoolean.True }))
  def greaterThan(rhs: WomValue): Try[WomBoolean] = invalid(s"$this > $rhs")
  def greaterThanOrEqual(rhs: WomValue): Try[WomBoolean] =
    Try(WomBoolean(Seq(greaterThan _, equals _).exists{ p => p(rhs).get == WomBoolean.True }))
  def or(rhs: WomValue): Try[WomBoolean] = invalid(s"$this || $rhs")
  def and(rhs: WomValue): Try[WomBoolean] = invalid(s"$this && $rhs")
  def not: Try[WomValue] = invalid(s"!$this")
  def unaryPlus: Try[WomValue] = invalid(s"+$this")
  def unaryMinus: Try[WomValue] = invalid(s"-$this")
  def typeName: String = womType.getClass.getSimpleName

  def toWomString: String = throw new UnsupportedOperationException(s"$getClass does not implement toWomString")

  /* This emits the value as a string.  In other words, the String value that
   * would be inserted into the command line.
   *
   * WomString("foobar") -> foobar
   *
   * toWomString is a good approximate implementation, though not sufficient
   * for types like WomString where extra syntax is added on
   */
  def valueString: String = toWomString

  def collectAsSeq[T <: WomValue](filterFn: PartialFunction[WomValue, T]): Seq[T] = {
    if (filterFn.isDefinedAt(this)) Seq(filterFn(this)) else Nil
  }

  private def symbolHash(hash: String) = SymbolHash((this.getClass.getCanonicalName + hash).md5Sum)

  private def symbolHash[K](hashedMap: Map[K, SymbolHash])(implicit ord: Ordering[K]): SymbolHash = {
    // productIterator returns an Iterator over the elements of a Tuple2 Map entry.
    val concatenatedMap = TreeMap(hashedMap.toArray: _*) flatMap { _.productIterator } mkString ""
    symbolHash(concatenatedMap)
  }

  def computeHash(implicit hasher: FileHasher): SymbolHash = {
    this match {
      case w: WomObject => symbolHash(w.values safeMapValues { _.computeHash(hasher) })
      case w: WomMap => symbolHash(w.value map { case (k, v) => k.computeHash(hasher) -> v.computeHash(hasher) })
      case w: WomArray => symbolHash(w.value map { _.computeHash(hasher) } mkString "")
      case w: WomFile => hasher(w)
      case w => symbolHash(w.valueString)
    }
  }

  def asWomExpression: ValueAsAnExpression = ValueAsAnExpression(this)

  /**
    * Perform any potentially async initialization on this wom value before it can be used to evaluate an expression
    * or instantiate a command for instance.
    * 
    * TODO: It would be better if the return type was the concrete one instead of a generic WomValue, but this 
    * seems hard to do without WomValue being parameterized.
    */
  def initialize(ioFunctionSet: IoFunctionSet): IOChecked[WomValue] = IOChecked.pure(this)
}

object WomValue {
  /**
    * Returns the womValue with all collections recursively limited to maximum length `maxElements`.
    *
    * @param womValue    The original womValue.
    * @param maxElements The maximum number of elements per collection.
    * @return The womValue with maximum maxElements per collection.
    */
  def takeMaxElements(womValue: WomValue, maxElements: Int): WomValue = {
    def takeMaxElements(recursiveWomValue: WomValue): WomValue = {
      recursiveWomValue match {
        case WomArray(womType, values) =>
          val subset = values.take(maxElements)
          WomArray(womType, subset map takeMaxElements)
        case WomMap(womType, values) =>
          val subset = values.take(maxElements)
          WomMap(
            womType,
            subset map {
              case (mapKey, mapValue) => takeMaxElements(mapKey) -> takeMaxElements(mapValue)
            }
          )
        case objectLike: WomObjectLike =>
          // First take only a limited number of the top-level elements.
          val shallowSubset = objectLike.values.take(maxElements)
          // Then recursively take only a limited number of elements.
          val deepSubset = shallowSubset map {
            case (mapKey, mapValue) => mapKey -> takeMaxElements(mapValue)
          }
          objectLike.copyWith(deepSubset)
        case WomOptionalValue(innerType, valueOption) =>
          WomOptionalValue(innerType, valueOption map takeMaxElements)
        case WomPair(left, right) => WomPair(takeMaxElements(left), takeMaxElements(right))
        case _ => recursiveWomValue
      }
    }

    takeMaxElements(womValue)
  }
}
