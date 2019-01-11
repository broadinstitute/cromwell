package cromwell.services.womtool.models

import io.circe.{Decoder, Encoder}
import io.circe.generic.semiauto._
import wom.types.WomType
import cromwell.services.womtool.models.WomTypeJsonSupport.womTypeEncoder // IntelliJ is lying
import cromwell.services.womtool.models.WomTypeJsonSupport.womTypeDecoder

case class InputDescription(name: String, valueType: WomType, optional: Boolean)

object InputDescription {
  implicit val inputDescriptionEncoder: Encoder[InputDescription] = deriveEncoder[InputDescription]
  implicit val inputDescriptionDecoder: Decoder[InputDescription] = deriveDecoder[InputDescription]
}
