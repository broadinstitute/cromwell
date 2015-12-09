package cromwell.binding.values

import cromwell.binding.{SymbolHash, FileHasher}
import cromwell.binding.types.{WdlStringType, WdlArrayType, WdlObjectType, WdlPrimitiveType}
import cromwell.util.StringUtil._

import scala.language.postfixOps
import scala.util.{Failure, Success, Try}

object WdlArray {
  def fromTsv(tsv: String): WdlArray = {
    WdlArray(WdlArrayType(WdlArrayType(WdlStringType)), tsv.replaceAll("[\r\n]+$", "").split("[\n\r]").toSeq map { line =>
      WdlArray(WdlArrayType(WdlStringType), line.split("\t").toSeq.map(WdlString))
    })
  }
}

case class WdlArray(wdlType: WdlArrayType, value: Seq[WdlValue]) extends WdlValue {
  val typesUsedInValue = Set(value map {_.wdlType}: _*)
  if (typesUsedInValue.size == 1 && typesUsedInValue.head != wdlType.memberType) {
    throw new UnsupportedOperationException(s"Could not construct array of type $wdlType with this value: $value")
  }
  if (typesUsedInValue.size > 1) {
    throw new UnsupportedOperationException(s"Cannot construct array with a mixed types: $value")
  }

  override def toWdlString: String = s"[${value.map(_.toWdlString).mkString(", ")}]"
  override def toString = toWdlString

  def map[R <: WdlValue](f: WdlValue => R): WdlArray = {
    value.map{f} match {
      case s:Seq[R] if s.nonEmpty => WdlArray(WdlArrayType(s.head.wdlType), s)
      case _ => this
    }
  }

  def tsvSerialize: Try[String] = {
    wdlType.memberType match {
      case t: WdlPrimitiveType => Success(value.map(_.valueString).mkString("\n"))
      case WdlObjectType => WdlObject.tsvSerializeArray(value map { _.asInstanceOf[WdlObject] })
      case _ => Failure(new UnsupportedOperationException("Can only TSV serialize an Array[Primitive]"))
    }
  }

  override def collectAsSeq[T <: WdlValue](filterFn: PartialFunction[WdlValue, T]): Seq[T] = {
    value flatMap { _.collectAsSeq(filterFn) }
  }

  override def getHash(implicit hasher: FileHasher): SymbolHash = {
    val hashedArray = value map { _.getHash }
    val accumulatedHash = hashedArray mkString ""
    SymbolHash((getClass.getCanonicalName+accumulatedHash).md5Sum)
  }
}
