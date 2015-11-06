package cromwell.binding.values

import cromwell.binding.WdlExpressionException
import cromwell.binding.types.WdlType
import cromwell.engine.{FileHasher, Hash}
import cromwell.util.StringDigestion._

import scala.util.{Failure, Try}

trait WdlValue {
  val wdlType: WdlType
  def invalid(operation: String) = Failure(new WdlExpressionException(s"Cannot perform operation: $operation"))
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

  /**
   * @return MD5 Hash.
   */
  def getHash(implicit hasher: FileHasher): Hash = valueString.md5Sum
}
