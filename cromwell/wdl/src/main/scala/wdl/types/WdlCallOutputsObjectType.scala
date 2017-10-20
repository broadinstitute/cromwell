package wdl.types

import wdl.WdlCall
import wdl.values.WdlCallOutputsObject
import wom.types.WomType

case class WdlCallOutputsObjectType(call: WdlCall) extends WomType {
  val toDisplayString: String = "Object"

  override protected def coercion = {
    case o: WdlCallOutputsObject => o
  }
}
