package wom.values

import wom.types.WomEnumerationType

case class WomEnumerationValue(womType: WomEnumerationType, value: String) extends WomValue {

  override def toWomString: String = value
}
