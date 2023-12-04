package wom.values

import wom.types.WomPairType

case class WomPair(left: WomValue, right: WomValue) extends WomValue {
  override val womType = WomPairType(left.womType, right.womType)

  override def toWomString = s"(${left.toWomString}, ${right.toWomString})"

  override def collectAsSeq[T <: WomValue](filterFn: PartialFunction[WomValue, T]): Seq[T] =
    left.collectAsSeq(filterFn) ++ right.collectAsSeq(filterFn)
}
