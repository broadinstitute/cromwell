package wdl.types

import wdl.WdlCall
import wdl.values.WdlCallOutputsObject
import wom.types.WdlType

case class WdlCallOutputsObjectType(call: WdlCall) extends WdlType {
  val toWdlString: String = "Object"

  override protected def coercion = {
    case o: WdlCallOutputsObject => o
  }
}
