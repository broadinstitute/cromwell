package wdl4s.values

import lenthall.util.TryUtil
import wdl4s.TsvSerializable
import wdl4s.types._
import wdl4s.util.FileUtil
import wdl4s.values.WdlArray.WdlArrayLike

import scala.language.postfixOps
import scala.util.{Failure, Success, Try}

object WdlMap {
  def coerceMap(m: Map[_, _], wdlMapType: WdlMapType): WdlMap = {
    val coerced = m map { case(k, v) => wdlMapType.keyType.coerceRawValue(k) -> wdlMapType.valueType.coerceRawValue(v) }
    val failures = coerced flatMap { case(k,v) => Seq(k,v) } collect { case f:Failure[_] => f }
    failures match {
      case f: Iterable[Failure[_]] if f.nonEmpty =>
        throw new UnsupportedOperationException(s"Failed to coerce one or more keys or values for creating a ${wdlMapType.toWdlString}:\n${TryUtil.stringifyFailures(f)}}")
      case _ =>
        val mapCoerced = coerced map { case (k, v) => k.get -> v.get }

        val keyType = WdlType.homogeneousTypeFromValues(mapCoerced map { case (k, v) => k })
        val valueType = WdlType.homogeneousTypeFromValues(mapCoerced map { case (k, v) => v })

        WdlMap(WdlMapType(keyType, valueType), mapCoerced)
    }
  }

  def fromTsv(tsv: String, wdlMapType: WdlMapType = WdlMapType(WdlAnyType, WdlAnyType)): Try[WdlMap] = {
    FileUtil.parseTsv(tsv) match {
      case Success(table) if table.isEmpty => Success(WdlMap(wdlMapType, Map.empty[WdlValue, WdlValue]))
      case Success(table) if table.head.length != 2 => Failure(new UnsupportedOperationException("TSV must be 2 columns to convert to a Map"))
      case Success(table) => Try(coerceMap(table.map(row => row(0) -> row(1)).toMap, wdlMapType))
      case Failure(e) => Failure(e)
    }
  }

  def apply(m: Map[WdlValue, WdlValue]): WdlMap = {
    val keyType = WdlType.lowestCommonSubtype(m.keys.map(_.wdlType))
    val valueType = WdlType.lowestCommonSubtype(m.keys.map(_.wdlType))
    WdlMap(WdlMapType(keyType, valueType), m)
  }
}

case class WdlMap(wdlType: WdlMapType, value: Map[WdlValue, WdlValue]) extends WdlValue with WdlArrayLike with TsvSerializable {
  val typesUsedInKey = value.map { case (k,v) => k.wdlType }.toSet

  if (typesUsedInKey.size == 1 && typesUsedInKey.head != wdlType.keyType)
    throw new UnsupportedOperationException(s"Could not construct a $wdlType as this value: $value")

  if (typesUsedInKey.size > 1)
    throw new UnsupportedOperationException(s"Cannot construct $wdlType with mixed types: $value")

  val typesUsedInValue = value.map { case (k,v) => v.wdlType }.toSet

  if (typesUsedInValue.size == 1 && typesUsedInValue.head != wdlType.valueType)
    throw new UnsupportedOperationException(s"Could not construct a $wdlType as this value: $value")

  if (typesUsedInValue.size > 1)
    throw new UnsupportedOperationException(s"Cannot construct $wdlType with mixed types: $value")

  override def toWdlString: String =
    "{" + value.map {case (k,v) => s"${k.toWdlString}: ${v.toWdlString}"}.mkString(", ") + "}"

  def tsvSerialize: Try[String] = {
    (wdlType.keyType, wdlType.valueType) match {
      case (wdlTypeKey: WdlPrimitiveType, wdlTypeValue: WdlPrimitiveType) =>
        Success(value.map({case (k, v) => s"${k.valueString}\t${v.valueString}"}).mkString("\n"))
      case _ =>
        Failure(new UnsupportedOperationException("Can only TSV serialize a Map[Primitive, Primitive]"))
    }
  }

  def map(f: PartialFunction[((WdlValue, WdlValue)), (WdlValue, WdlValue)]): WdlMap = {
    value map f match {
      case m: Map[WdlValue, WdlValue] if m.nonEmpty => WdlMap(WdlMapType(m.head._1.wdlType, m.head._2.wdlType), m)
      case _ => this
    }
  }

  override def collectAsSeq[T <: WdlValue](filterFn: PartialFunction[WdlValue, T]): Seq[T] = {
    val collected = value flatMap {
      case (k, v) => Seq(k.collectAsSeq(filterFn), v.collectAsSeq(filterFn))
    }
    collected.flatten.toSeq
  }

  // For WdlArrayLike:
  override lazy val arrayType: WdlArrayType = WdlArrayType(WdlPairType(wdlType.keyType, wdlType.valueType))
  override lazy val asArray: WdlArray = WdlArray(arrayType, value.toSeq map { case (k, v) => WdlPair(k, v) })
}
