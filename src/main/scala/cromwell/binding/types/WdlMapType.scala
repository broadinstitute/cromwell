package cromwell.binding.types

import cromwell.binding.values._
import spray.json.JsObject

case class WdlMapType(keyType: WdlType, valueType: WdlType) extends WdlType {
  val toWdlString: String = s"Map[${keyType.toWdlString}, ${valueType.toWdlString}]"

  override protected def coercion = {
    case m: Map[_, _] if m.nonEmpty => WdlMap.coerceMap(m, this)
    case m: Map[_, _] if m.isEmpty => WdlMap(WdlMapType(keyType, valueType), Map())
    case js: JsObject if js.fields.nonEmpty => WdlMap.coerceMap(js.fields, this)
    case wdlMap: WdlMap => WdlMap.coerceMap(wdlMap.value, this)
  }
}
