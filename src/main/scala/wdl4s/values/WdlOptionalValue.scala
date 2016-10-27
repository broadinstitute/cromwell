package wdl4s.values

import wdl4s.types.{WdlOptionalType, WdlType}

case class WdlOptionalValue(innerType: WdlType, value: Option[WdlValue]) extends WdlValue {
  override val wdlType = WdlOptionalType(innerType)
  override val toWdlString = value map { _.toWdlString } getOrElse "null"
}

object WdlOptionalValue {

  def apply(value: WdlValue): WdlOptionalValue = Option(value) match {
    case someValue @ Some(innerValue) => WdlOptionalValue(innerValue.wdlType, someValue)
    case None => throw new NullPointerException(s"Cannot use apply for a null WdlValue")
  }

  def none(wdlType: WdlType) = WdlOptionalValue(wdlType, None)
}
