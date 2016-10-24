package wdl4s.values

import wdl4s.types.WdlPairType

case class WdlPair(left: WdlValue, right: WdlValue) extends WdlValue {
  override val wdlType = WdlPairType(left.wdlType, right.wdlType)

  override def toWdlString = s"(${left.toWdlString}, ${right.toWdlString})"
}
