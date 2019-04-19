package wom.values

import cats.Applicative
import cats.syntax.functor._
import common.util.TryUtil
import wom.TsvSerializable
import wom.types._
import wom.util.FileUtil
import wom.values.WomArray.WomArrayLike
import cats.syntax.traverse._
import cats.instances.list._
import common.validation.IOChecked.IOChecked
import wom.expression.IoFunctionSet

import scala.language.higherKinds
import scala.util.{Failure, Success, Try}

object WomMap {
  def coerceMap(m: Map[_, _], womMapType: WomMapType): WomMap = {
    val coerced: Map[Try[WomValue], Try[WomValue]] = m map { case(k, v) => womMapType.keyType.coerceRawValue(k) -> womMapType.valueType.coerceRawValue(v) }
    val failures = coerced flatMap { case(k,v) => Seq(k,v) } collect { case f:Failure[_] => f }
    failures match {
      case f: Iterable[Failure[_]] if f.nonEmpty =>
        throw new UnsupportedOperationException(s"Failed to coerce one or more keys or values for creating a ${womMapType.stableName}:\n${TryUtil.stringifyFailures(f)}}")
      case _ =>
        val mapCoerced = coerced map { case (k, v) => k.get -> v.get }

        val keyType = WomType.homogeneousTypeFromValues(mapCoerced map { case (k, _) => k })
        val valueType = WomType.homogeneousTypeFromValues(mapCoerced map { case (_, v) => v })

        WomMap(WomMapType(keyType, valueType), mapCoerced)
    }
  }

  def fromTsv(tsv: String, womMapType: WomMapType = WomMapType(WomAnyType, WomAnyType)): Try[WomMap] = {
    FileUtil.parseTsv(tsv) match {
      case Success(table) if table.isEmpty => Success(WomMap(womMapType, Map.empty[WomValue, WomValue]))
      case Success(table) if table.head.length != 2 => Failure(new UnsupportedOperationException("TSV must be 2 columns to convert to a Map"))
      case Success(table) => Try(coerceMap(table.map(row => row(0) -> row(1)).toMap, womMapType))
      case Failure(e) => Failure(e)
    }
  }

  def apply(m: Map[WomValue, WomValue]): WomMap = {
    val keyType = WomType.lowestCommonSubtype(m.keys.map(_.womType))
    val valueType = WomType.lowestCommonSubtype(m.values.map(_.womType))
    coerceMap(m, WomMapType(keyType, valueType))
  }
}

final case class WomMap private(womType: WomMapType, value: Map[WomValue, WomValue]) extends WomValue with WomArrayLike with TsvSerializable {
  val typesUsedInKey = value.map { case (k, _) => k.womType }.toSet

  if (typesUsedInKey.size == 1 && typesUsedInKey.head != womType.keyType)
    throw new UnsupportedOperationException(s"Could not construct a $womType with map keys of unexpected type: [${value.keys.mkString(", ")}]")

  if (typesUsedInKey.size > 1)
    throw new UnsupportedOperationException(s"Cannot construct $womType with mixed types in map keys: [${value.keys.mkString(", ")}]")

  val typesUsedInValue = value.map { case (_, v) => v.womType }.toSet

  if (typesUsedInValue.size == 1 && typesUsedInValue.head != womType.valueType)
    throw new UnsupportedOperationException(s"Could not construct a $womType with map values of unexpected type: [${value.values.mkString(", ")}]")

  if (typesUsedInValue.size > 1)
    throw new UnsupportedOperationException(s"Cannot construct $womType with mixed types in map values: [${value.values.mkString(", ")}]")

  override def toWomString: String =
    "{" + value.map {case (k,v) => s"${k.toWomString}: ${v.toWomString}"}.mkString(", ") + "}"

  def tsvSerialize: Try[String] = {
    (womType.keyType, womType.valueType) match {
      case (_: WomPrimitiveType, _: WomPrimitiveType) =>
        Success(value.map({case (k, v) => s"${k.valueString}\t${v.valueString}"}).mkString("\n"))
      case _ =>
        Failure(new UnsupportedOperationException("Can only TSV serialize a Map[Primitive, Primitive]"))
    }
  }

  def map(f: PartialFunction[(WomValue, WomValue), (WomValue, WomValue)]): WomMap = {
    value map f match {
      case m: Map[WomValue, WomValue] if m.nonEmpty => WomMap(WomMapType(m.head._1.womType, m.head._2.womType), m)
      case _ => this
    }
  }

  def traverseValues[R <: WomValue, G[_]](f: WomValue => G[R])(implicit applicative: Applicative[G]): G[WomMap] = {
    if (value.isEmpty) applicative.pure(this)
    else {
      val traverseFunction: (WomValue, WomValue) => G[(WomValue, R)] = {
        case (key, v) => applicative.map(f(v)) { key -> _ }
      }
      applicative.map(value.toList.traverse[G, (WomValue, R)](traverseFunction.tupled)) { mapped =>
        WomMap(mapped.toMap)
      }
    }
  }

  override def collectAsSeq[T <: WomValue](filterFn: PartialFunction[WomValue, T]): Seq[T] = {
    val collected = value flatMap {
      case (k, v) => Seq(k.collectAsSeq(filterFn), v.collectAsSeq(filterFn))
    }
    collected.flatten.toSeq
  }

  override def initialize(ioFunctionSet: IoFunctionSet): IOChecked[WomValue] = traverseValues(_.initialize(ioFunctionSet)).widen

  // For WomArrayLike:
  override lazy val arrayType: WomArrayType = WomArrayType(WomPairType(womType.keyType, womType.valueType))
  override lazy val asArray: WomArray = WomArray(arrayType, value.toSeq map { case (k, v) => WomPair(k, v) })
}
