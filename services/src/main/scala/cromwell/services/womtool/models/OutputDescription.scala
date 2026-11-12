package cromwell.services.womtool.models

import cromwell.services.womtool.models.WomTypeJsonSupport.womTypeEncoder
import spray.json.{JsObject, JsonWriter, JsString, JsValue}
import wom.types.WomType

case class OutputDescription(name: String, valueType: WomType, typeDisplayName: String)

object OutputDescription {
  implicit val outputDescriptionWriter: JsonWriter[OutputDescription] = new JsonWriter[OutputDescription] {
    override def write(a: OutputDescription): JsValue = JsObject(
      "name" -> JsString(a.name),
      "valueType" -> womTypeEncoder.write(a.valueType),
      "typeDisplayName" -> JsString(a.typeDisplayName)
    )
  }
}
