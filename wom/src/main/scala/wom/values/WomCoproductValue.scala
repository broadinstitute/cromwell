package wom.values

import wom.types.{WomCoproductType, WomType}

case class WomCoproductValue(innerType: WomCoproductType, womValue: WomValue) extends WomValue {
  override def womType: WomType = innerType
}
