package cromwell.binding.values

import cromwell.binding.types.WdlArrayType

case class WdlArray(wdlType: WdlArrayType, value: Seq[WdlValue]) extends WdlValue {
  val typesUsedInValue = Set(value map {_.wdlType}: _*)
  if (typesUsedInValue.size == 1 && typesUsedInValue.head != wdlType.memberType) {
    throw new UnsupportedOperationException(s"Could not construct array of type $wdlType with this value: $value")
  }
  if (typesUsedInValue.size > 1) {
    throw new UnsupportedOperationException(s"Cannot construct array with a mixed types: $value")
  }

  override def toWdlString: String = s"[${value.map(_.toWdlString).mkString(", ")}]"

  def map[R <: WdlValue](f: WdlValue => R): WdlArray = {
    value.map{f} match {
      case s:Seq[R] if s.nonEmpty => WdlArray(WdlArrayType(s.head.wdlType), s)
      case _ => this
    }
  }
}