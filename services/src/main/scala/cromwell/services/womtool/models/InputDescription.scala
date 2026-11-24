package cromwell.services.womtool.models

import cromwell.services.womtool.models.WomExpressionJsonSupport.womExpressionEncoder
import cromwell.services.womtool.models.WomTypeJsonSupport.womTypeEncoder
import spray.json.{DefaultJsonProtocol, JsonFormat, JsonWriter, NullOptions, RootJsonFormat}
import wom.expression.WomExpression
import wom.types.WomType

case class InputDescription(name: String,
                            valueType: WomType,
                            typeDisplayName: String,
                            optional: Boolean,
                            default: Option[WomExpression]
)

object InputDescription {
  // NullOptions writes a missing `default` as `null`, which the describe endpoint has always returned.
  private object Protocol extends DefaultJsonProtocol with NullOptions {
    implicit val womTypeFormat: JsonFormat[WomType] = lift(womTypeEncoder)
    implicit val womExpressionFormat: JsonFormat[WomExpression] = lift(womExpressionEncoder)
    val inputDescriptionFormat: RootJsonFormat[InputDescription] = jsonFormat5(InputDescription.apply)
  }

  implicit val inputDescriptionWriter: JsonWriter[InputDescription] = Protocol.inputDescriptionFormat
}
