package wom.values

import wom.types.WomCoproductType

case class WomCoproductValue(womType: WomCoproductType, womValue: WomValue) extends WomValue {

  override def toWomString: String = {
    womValue.toWomString
  }
}
