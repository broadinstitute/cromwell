package cromwell.services.womtool.models

import cromwell.services.womtool.models.WomTypeJsonSupport.womTypeEncoder // IntelliJ is lying
import cromwell.services.womtool.models.WomExpressionJsonSupport.womExpressionEncoder
import io.circe.generic.JsonCodec
import wom.expression.WomExpression
import wom.types.WomType

@JsonCodec(encodeOnly=true)
case class InputDescription(name: String, valueType: WomType, typeDisplayName: String, optional: Boolean, default: Option[WomExpression])
