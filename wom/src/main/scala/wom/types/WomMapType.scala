package wom.types

import spray.json.JsObject
import wom.values.{WomMap, WomObjectLike}

case class WomMapType(keyType: WomType, valueType: WomType) extends WomType {
  val stableName: String = s"Map[${keyType.stableName}, ${valueType.stableName}]"

  override protected def coercion = {
    case m: Map[_, _] if m.nonEmpty => WomMap.coerceMap(m, this)
    case m: Map[_, _] if m.isEmpty => WomMap(WomMapType(keyType, valueType), Map())
    case js: JsObject if js.fields.nonEmpty => WomMap.coerceMap(js.fields, this)
    case womMap: WomMap => WomMap.coerceMap(womMap.value, this)
    case o: WomObjectLike => WomMap.coerceMap(o.values, this)
  }

  override def typeSpecificIsCoerceableFrom(otherType: WomType): Boolean = otherType match {
    case m: WomMapType => keyType.isCoerceableFrom(m.keyType) && valueType.isCoerceableFrom(m.valueType)
    case _: WomObjectTypeLike => keyType.isCoerceableFrom(WomStringType) && valueType.isCoerceableFrom(WomStringType)
    case _ => false
  }

  def equivalentArrayType: WomArrayType = WomArrayType(WomPairType(keyType, valueType))
}
