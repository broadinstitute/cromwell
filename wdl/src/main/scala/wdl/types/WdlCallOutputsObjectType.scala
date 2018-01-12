package wdl.types

import wdl.WdlCall
import wdl.values.WdlCallOutputsObject
import wom.types.WomObjectTypeLike

case class WdlCallOutputsObjectType(call: WdlCall) extends WomObjectTypeLike {
  val toDisplayString: String = "Object"

  override protected def coercion = {
    case o: WdlCallOutputsObject => o
  }
}
