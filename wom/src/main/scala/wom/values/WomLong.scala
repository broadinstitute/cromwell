package wom.values

import wom.types.WomLongType

case class WomLong(value: Long) extends WomPrimitive {
  val womType = WomLongType

  override def toWomString = value.toString
}
