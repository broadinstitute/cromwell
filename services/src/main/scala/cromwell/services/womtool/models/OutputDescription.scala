package cromwell.services.womtool.models

import cromwell.services.womtool.models.WomTypeJsonSupport.womTypeEncoder
import spray.json.{DefaultJsonProtocol, JsonFormat, JsonWriter, RootJsonFormat}
import wom.types.WomType

case class OutputDescription(name: String, valueType: WomType, typeDisplayName: String)

object OutputDescription {
  private object Protocol extends DefaultJsonProtocol {
    implicit val womTypeFormat: JsonFormat[WomType] = lift(womTypeEncoder)
    val outputDescriptionFormat: RootJsonFormat[OutputDescription] = jsonFormat3(OutputDescription.apply)
  }

  implicit val outputDescriptionWriter: JsonWriter[OutputDescription] = Protocol.outputDescriptionFormat
}
