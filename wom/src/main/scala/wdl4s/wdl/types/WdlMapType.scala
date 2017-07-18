package wdl4s.wdl.types

import wdl4s.wdl.values._
import spray.json.JsObject

case class WdlMapType(keyType: WdlType, valueType: WdlType) extends WdlType {
  val toWdlString: String = s"Map[${keyType.toWdlString}, ${valueType.toWdlString}]"

  override protected def coercion = {
    case m: Map[_, _] if m.nonEmpty => WdlMap.coerceMap(m, this)
    case m: Map[_, _] if m.isEmpty => WdlMap(WdlMapType(keyType, valueType), Map())
    case js: JsObject if js.fields.nonEmpty => WdlMap.coerceMap(js.fields, this)
    case wdlMap: WdlMap => WdlMap.coerceMap(wdlMap.value, this)
    case o: WdlObject => WdlMap.coerceMap(o.value, this)
  }

  override def isCoerceableFrom(otherType: WdlType): Boolean = otherType match {
    case m: WdlMapType => keyType.isCoerceableFrom(m.keyType) && valueType.isCoerceableFrom(m.valueType)
    case WdlObjectType => keyType.isCoerceableFrom(WdlStringType) && valueType.isCoerceableFrom(WdlStringType)
    case _ => false
  }
}
