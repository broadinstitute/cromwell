package wdl.draft2.types

import wdl.draft2.WdlCall
import wdl.draft2.values.WdlCallOutputsObject
import wom.types.WomObjectTypeLike

case class WdlCallOutputsObjectType(call: WdlCall) extends WomObjectTypeLike {
  val toDisplayString: String = "Object"

  override protected def coercion = {
    case o: WdlCallOutputsObject => o
  }
}
