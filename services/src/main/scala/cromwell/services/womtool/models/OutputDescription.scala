package cromwell.services.womtool.models

import io.circe.Encoder
import io.circe.generic.semiauto.deriveEncoder
import wom.types.WomType
import cromwell.services.womtool.models.WomTypeJsonSupport.womTypeEncoder // IntelliJ is lying

case class OutputDescription(name: String, valueType: WomType, typeDisplayName: String)

object OutputDescription {
  implicit val outputDescriptionEncoder: Encoder[OutputDescription] = deriveEncoder[OutputDescription]
}
