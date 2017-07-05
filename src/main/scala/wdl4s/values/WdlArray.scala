package wdl4s.values

import lenthall.util.TryUtil
import wdl4s.TsvSerializable
import wdl4s.types._
import wdl4s.values.WdlArray.WdlArrayLike

import scala.language.postfixOps
import scala.util.{Failure, Success, Try}

object WdlArray {
  def fromTsv(tsv: String): WdlArray = {
    WdlArray(WdlArrayType(WdlArrayType(WdlStringType)), tsv.replaceAll("[\r\n]+$", "").split("[\n\r]+").toSeq map { line =>
      WdlArray(WdlArrayType(WdlStringType), line.split("\t").toSeq.map(WdlString))
    })
  }

  def apply(wdlType: WdlArrayType, value: Seq[WdlValue]): WdlArray = {
    if (wdlType == WdlMaybeEmptyArrayType.EmptyArrayType && value.nonEmpty) {
      throw new UnsupportedOperationException(s"An ${wdlType.toWdlString} must be empty but instead has value: ${value.mkString(", ")}")
    }
    if (wdlType.guaranteedNonEmpty && value.isEmpty) {
      throw new UnsupportedOperationException(s"An ${wdlType.toWdlString} must contain at least one element")
    }

    val coercedValue = TryUtil.sequence(value map wdlType.memberType.coerceRawValue)
    coercedValue match {
      case Success(coercedArray) => new WdlArray(wdlType, coercedArray) {}
      case Failure(f) => throw new UnsupportedOperationException(s"Could not construct array of type $wdlType with this value: $value", f)
    }
  }

  trait WdlArrayLike {
    def arrayType: WdlArrayType
    def asArray: WdlArray
  }

  object WdlArrayLike {
    def unapply(v: WdlValue): Option[WdlArray] = v match {
      case wal: WdlArrayLike => Option(wal.asArray)
      case _ => None
    }
  }
}

sealed abstract case class WdlArray(wdlType: WdlArrayType, value: Seq[WdlValue]) extends WdlValue with WdlArrayLike with TsvSerializable {

  val nonEmpty = value.nonEmpty
  override def toWdlString: String = s"[${value.map(_.toWdlString).mkString(", ")}]"
  override def toString = toWdlString

  def map[R <: WdlValue](f: WdlValue => R): WdlArray = {
    value.map{f} match {
      case s: Seq[R] if s.nonEmpty => WdlArray(WdlArrayType(s.head.wdlType), s)
      case _ => this
    }
  }

  def tsvSerialize: Try[String] = {
    wdlType.memberType match {
      case t: WdlPrimitiveType => Success(value.map(_.valueString).mkString(start = "", sep = "\n", end = "\n"))
      case WdlObjectType => WdlObject.tsvSerializeArray(value map { _.asInstanceOf[WdlObject] })
      case WdlArrayType(t: WdlPrimitiveType) =>
        val tsvString = value.collect({ case a: WdlArray => a }) map { a =>
          a.value.collect({ case p: WdlPrimitive => p.valueString }).mkString(start = "", sep = "\t", end = "\n")
        } mkString

        Success(tsvString)
      case _ => Failure(new UnsupportedOperationException(s"Cannot TSV serialize a ${this.wdlType.toWdlString} (valid types are Array[Primitive], Array[Array[Primitive]], or Array[Object])"))
    }
  }

  override def collectAsSeq[T <: WdlValue](filterFn: PartialFunction[WdlValue, T]): Seq[T] = {
    value flatMap { _.collectAsSeq(filterFn) }
  }

  // For WdlArrayLike:
  override val arrayType: WdlArrayType = wdlType
  override val asArray: WdlArray = this
}
