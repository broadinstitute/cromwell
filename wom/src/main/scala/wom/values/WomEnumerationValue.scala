package wom.values

import wom.types.WomEnumerationType

case class WomEnumerationValue(womType: WomEnumerationType, value: String) extends WomPrimitive {

  override def toWomString: String = value
}
