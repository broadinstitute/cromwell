package wom.types

import spray.json.JsObject
import wom.values.{WomMap, WomObject}

case class WomMapType(keyType: WomType, valueType: WomType) extends WomType {
  val toDisplayString: String = s"Map[${keyType.toDisplayString}, ${valueType.toDisplayString}]"

  override protected def coercion = {
    case m: Map[_, _] if m.nonEmpty => WomMap.coerceMap(m, this)
    case m: Map[_, _] if m.isEmpty => WomMap(WomMapType(keyType, valueType), Map())
    case js: JsObject if js.fields.nonEmpty => WomMap.coerceMap(js.fields, this)
    case womMap: WomMap => WomMap.coerceMap(womMap.value, this)
    case o: WomObject => WomMap.coerceMap(o.value, this)
  }

  override def isCoerceableFrom(otherType: WomType): Boolean = otherType match {
    case m: WomMapType => keyType.isCoerceableFrom(m.keyType) && valueType.isCoerceableFrom(m.valueType)
    case WomObjectType => keyType.isCoerceableFrom(WomStringType) && valueType.isCoerceableFrom(WomStringType)
    case _ => false
  }

  def equivalentArrayType: WomArrayType = WomArrayType(WomPairType(keyType, valueType))
}
