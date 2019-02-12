package wom.values

import cats.Applicative
import cats.instances.list._
import cats.syntax.traverse._
import cats.syntax.functor._
import common.util.TryUtil
import common.validation.IOChecked.IOChecked
import wom.TsvSerializable
import wom.expression.IoFunctionSet
import wom.types._
import wom.values.WomArray.WomArrayLike

import scala.language.{higherKinds, postfixOps}
import scala.util.{Failure, Success, Try}

object WomArray {
  def empty = WomArray(WomMaybeEmptyArrayType.EmptyArrayType, List.empty)

  def fromTsv(tsv: String): WomArray = {
    WomArray(WomArrayType(WomArrayType(WomStringType)), tsv.replaceAll("[\r\n]+$", "").split("[\n\r]+").toSeq map { line =>
      WomArray(WomArrayType(WomStringType), line.split("\t").toSeq.map(WomString))
    })
  }

  def apply(womType: WomArrayType, value: Seq[WomValue]): WomArray = {
    if (womType == WomMaybeEmptyArrayType.EmptyArrayType && value.nonEmpty) {
      throw new UnsupportedOperationException(s"An ${womType.stableName} must be empty but instead has value: ${value.mkString(", ")}")
    }
    if (womType.guaranteedNonEmpty && value.isEmpty) {
      throw new UnsupportedOperationException(s"An ${womType.stableName} must contain at least one element")
    }

    val coercedValue = TryUtil.sequence(value map womType.memberType.coerceRawValue)
    coercedValue match {
      case Success(coercedArray) => new WomArray(womType, coercedArray) {}
      case Failure(f) => throw new UnsupportedOperationException(s"Could not construct array of type $womType with this value: $value", f)
    }
  }

  def apply(value: Seq[WomValue]): WomArray = WomArray.apply(WomArrayType(WomType.homogeneousTypeFromValues(value)), value)

  trait WomArrayLike {
    def arrayType: WomArrayType
    def asArray: WomArray
  }

  object WomArrayLike {
    def unapply(v: WomValue): Option[WomArray] = v match {
      case wal: WomArrayLike => Option(wal.asArray)
      case _ => None
    }
  }
}

sealed abstract case class WomArray(womType: WomArrayType, value: Seq[WomValue]) extends WomValue with WomArrayLike with TsvSerializable {

  val nonEmpty = value.nonEmpty
  override def toWomString: String = s"[${value.map(_.toWomString).mkString(", ")}]"
  override def toString = toWomString

  def map[R <: WomValue](f: WomValue => R): WomArray = {
    value.map{f} match {
      case s: Seq[R] if s.nonEmpty => WomArray(WomArrayType(s.head.womType), s)
      case _ => this
    }
  }

  def traverse[R <: WomValue, G[_]](f: WomValue => G[R])(implicit applicative: Applicative[G]): G[WomArray] = {
    if (value.isEmpty) applicative.pure(this)
    else {
      applicative.map(value.toList.traverse(f)) { mapped =>
        WomArray(WomArrayType(mapped.head.womType), mapped)
      }
    }
  }

  override def initialize(ioFunctionSet: IoFunctionSet): IOChecked[WomValue] = traverse(_.initialize(ioFunctionSet)).widen

  def size = value.size

  def tsvSerialize: Try[String] = {
    womType.memberType match {
      case _: WomPrimitiveType => Success(value.map(_.valueString).mkString(start = "", sep = "\n", end = "\n"))
      case WomObjectType => WomObject.tsvSerializeArray(value map { _.asInstanceOf[WomObject] })
      case WomArrayType(_: WomPrimitiveType) =>
        val tsvString = value.collect({ case a: WomArray => a }) map { a =>
          a.value.collect({ case p: WomPrimitive => p.valueString }).mkString(start = "", sep = "\t", end = "\n")
        } mkString

        Success(tsvString)
      case _ => Failure(new UnsupportedOperationException(s"Cannot TSV serialize a ${this.womType.stableName} (valid types are Array[Primitive], Array[Array[Primitive]], or Array[Object])"))
    }
  }

  override def collectAsSeq[T <: WomValue](filterFn: PartialFunction[WomValue, T]): Seq[T] = {
    value flatMap { _.collectAsSeq(filterFn) }
  }

  // For WomArrayLike:
  override val arrayType: WomArrayType = womType
  override val asArray: WomArray = this
}
