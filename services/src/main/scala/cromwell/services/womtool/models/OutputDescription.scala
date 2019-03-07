package cromwell.services.womtool.models

import cromwell.services.womtool.models.WomTypeJsonSupport.womTypeEncoder // IntelliJ is lying
import io.circe.generic.JsonCodec
import wom.types.WomType

@JsonCodec(encodeOnly=true)
case class OutputDescription(name: String, valueType: WomType, typeDisplayName: String)
