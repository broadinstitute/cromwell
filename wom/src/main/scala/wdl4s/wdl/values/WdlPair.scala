package wdl4s.wdl.values

import wdl4s.wdl.types.WdlPairType

case class WdlPair(left: WdlValue, right: WdlValue) extends WdlValue {
  override val wdlType = WdlPairType(left.wdlType, right.wdlType)

  override def toWdlString = s"(${left.toWdlString}, ${right.toWdlString})"

  override def collectAsSeq[T <: WdlValue](filterFn: PartialFunction[WdlValue, T]): Seq[T] = {
    left.collectAsSeq(filterFn) ++ right.collectAsSeq(filterFn)
  }
}
