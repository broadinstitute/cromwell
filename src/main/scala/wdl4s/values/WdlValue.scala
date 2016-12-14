package wdl4s.values

import wdl4s.WdlExpressionException
import wdl4s.exception.VariableNotFoundException
import wdl4s.types.WdlType

import scala.collection.immutable.TreeMap
import scala.util.{Failure, Try}


trait WdlValue {
  val wdlType: WdlType
  def invalid(operation: String) = Failure(new WdlExpressionException(s"Cannot perform operation: $operation"))
  def emptyValue(value: WdlOptionalValue) = Failure(VariableNotFoundException(value.toString))
  def evaluateIfDefined[A <: WdlValue](optionalValue: WdlOptionalValue, operation: WdlValue => Try[A]): Try[A] = optionalValue match {
    case WdlOptionalValue(_, Some(v)) => operation(v)
    case _ => emptyValue(optionalValue)
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
  def toWdlString: String = ???

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
