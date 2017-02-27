package wdl4s.values

import wdl4s.WdlExpressionException
import wdl4s.exception.OptionalNotSuppliedException
import wdl4s.types.WdlType

import scala.collection.immutable.TreeMap
import scala.util.{Failure, Try}

trait WdlValue {
  val wdlType: WdlType
  def invalid(operation: String) = Failure(new WdlExpressionException(s"Cannot perform operation: $operation"))
  def emptyValueFailure(operationName: String) = Failure(OptionalNotSuppliedException(operationName))
  def evaluateIfDefined[A <: WdlValue](operationName: String, optionalValue: WdlOptionalValue, operation: WdlValue => Try[A]): Try[A] = optionalValue match {
    case WdlOptionalValue(_, Some(v)) => operation(v)
    case _ => emptyValueFailure(operationName)
  }
  
  def add(rhs: WdlValue): Try[WdlValue] = invalid(s"$this + $rhs")
  def subtract(rhs: WdlValue): Try[WdlValue] = invalid(s"$this - $rhs")
  def multiply(rhs: WdlValue): Try[WdlValue] = invalid(s"$this * $rhs")
  def divide(rhs: WdlValue): Try[WdlValue] = invalid(s"$this / $rhs")
  def mod(rhs: WdlValue): Try[WdlValue] = invalid(s"$this % $rhs")
  def equals(rhs: WdlValue): Try[WdlBoolean] = invalid(s"$this == $rhs")
  def notEquals(rhs: WdlValue): Try[WdlBoolean] = equals(rhs).map{x => WdlBoolean(!x.value)}
  def lessThan(rhs: WdlValue): Try[WdlBoolean] = invalid(s"$this < $rhs")
  def lessThanOrEqual(rhs: WdlValue): Try[WdlBoolean] =
    Try(WdlBoolean(Seq(lessThan _, equals _).exists{ p => p(rhs).get == WdlBoolean.True }))
  def greaterThan(rhs: WdlValue): Try[WdlBoolean] = invalid(s"$this > $rhs")
  def greaterThanOrEqual(rhs: WdlValue): Try[WdlBoolean] =
    Try(WdlBoolean(Seq(greaterThan _, equals _).exists{ p => p(rhs).get == WdlBoolean.True }))
  def or(rhs: WdlValue): Try[WdlBoolean] = invalid(s"$this || $rhs")
  def and(rhs: WdlValue): Try[WdlBoolean] = invalid(s"$this && $rhs")
  def not: Try[WdlValue] = invalid(s"!$this")
  def unaryPlus: Try[WdlValue] = invalid(s"+$this")
  def unaryMinus: Try[WdlValue] = invalid(s"-$this")
  def typeName: String = wdlType.getClass.getSimpleName

  /* This emits valid WDL source.  WdlString("foobar") -> "foobar" (quotes included) */
  def toWdlString: String = throw new NotImplementedError(s"$getClass does not implement toWdlString")

  /* This emits the value as a string.  In other words, the String value that
   * would be inserted into the command line.
   *
   * WdlString("foobar") -> foobar
   *
   * toWdlString is a good approximate implementation, though not sufficient
   * for types like WdlString where extra syntax is added on
   */
  def valueString: String = toWdlString

  def collectAsSeq[T <: WdlValue](filterFn: PartialFunction[WdlValue, T]): Seq[T] = {
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
      case w: WdlObject => symbolHash(w.value mapValues { _.computeHash(hasher) })
      case w: WdlMap => symbolHash(w.value map { case (k, v) => k.computeHash(hasher) -> v.computeHash(hasher) })
      case w: WdlArray => symbolHash(w.value map { _.computeHash(hasher) } mkString "")
      case w: WdlFile => hasher(w)
      case w => symbolHash(w.valueString)
    }
  }
}

object WdlValue {
  /**
    * Returns the wdlValue with all collections recursively limited to maximum length `maxElements`.
    *
    * @param wdlValue    The original wdlValue.
    * @param maxElements The maximum number of elements per collection.
    * @return The wdlValue with maximum maxElements per collection.
    */
  def takeMaxElements(wdlValue: WdlValue, maxElements: Int): WdlValue = {
    def takeMaxElements(recursiveWdlValue: WdlValue): WdlValue = {
      recursiveWdlValue match {
        case WdlArray(wdlType, values) =>
          val subset = values.take(maxElements)
          WdlArray(wdlType, subset map takeMaxElements)
        case WdlMap(wdlType, values) =>
          val subset = values.take(maxElements)
          WdlMap(
            wdlType,
            subset map {
              case (mapKey, mapValue) => takeMaxElements(mapKey) -> takeMaxElements(mapValue)
            }
          )
        case WdlObject(values) =>
          val subset = values.take(maxElements)
          WdlObject(subset map {
            case (mapKey, mapValue) => mapKey -> takeMaxElements(mapValue)
          })
        case WdlCallOutputsObject(call, values) =>
          val subset = values.take(maxElements)
          WdlCallOutputsObject(call, subset map {
            case (mapKey, mapValue) => mapKey -> takeMaxElements(mapValue)
          })
        case WdlOptionalValue(innerType, valueOption) =>
          WdlOptionalValue(innerType, valueOption map takeMaxElements)
        case WdlPair(left, right) => WdlPair(takeMaxElements(left), takeMaxElements(right))
        case _ => recursiveWdlValue
      }
    }

    takeMaxElements(wdlValue)
  }
}
