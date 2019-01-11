package cromwell.services.womtool.models

import wom.types.WomType

case class InputDescription(name: String, valueType: WomType, optional: Boolean)