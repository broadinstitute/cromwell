package cromwell.services.womtool.models

import cromwell.services.womtool.models.WomExpressionJsonSupport.womExpressionEncoder
import cromwell.services.womtool.models.WomTypeJsonSupport.womTypeEncoder
import spray.json.{JsBoolean, JsNull, JsObject, JsString, JsValue, JsonWriter}
import wom.expression.WomExpression
import wom.types.WomType

case class InputDescription(name: String,
                            valueType: WomType,
                            typeDisplayName: String,
                            optional: Boolean,
                            default: Option[WomExpression]
)

object InputDescription {
  implicit val inputDescriptionWriter: JsonWriter[InputDescription] = new JsonWriter[InputDescription] {
    override def write(a: InputDescription): JsValue = JsObject(
      "name" -> JsString(a.name),
      "valueType" -> womTypeEncoder.write(a.valueType),
      "typeDisplayName" -> JsString(a.typeDisplayName),
      "optional" -> JsBoolean(a.optional),
      "default" -> a.default.map(womExpressionEncoder.write).getOrElse(JsNull)
    )
  }
}
