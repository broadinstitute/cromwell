package cromwell.services.womtool.models

import io.circe.{Decoder, Encoder}
import io.circe.generic.semiauto.{deriveDecoder, deriveEncoder}
import wom.types.WomType
import cromwell.services.womtool.models.WomTypeJsonSupport.womTypeEncoder // IntelliJ is lying
import cromwell.services.womtool.models.WomTypeJsonSupport.womTypeDecoder

case class OutputDescription(name: String, valueType: WomType, typeDisplayName: String)

object OutputDescription {
  implicit val outputDescriptionEncoder: Encoder[OutputDescription] = deriveEncoder[OutputDescription]
  implicit val outputDescriptionDecoder: Decoder[OutputDescription] = deriveDecoder[OutputDescription]
}
