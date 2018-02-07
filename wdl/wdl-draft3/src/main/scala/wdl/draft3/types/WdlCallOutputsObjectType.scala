package wdl.draft3.types

import wdl.draft3.WdlCall
import wdl.draft3.values.WdlCallOutputsObject
import wom.types.WomObjectTypeLike

case class WdlCallOutputsObjectType(call: WdlCall) extends WomObjectTypeLike {
  val toDisplayString: String = "Object"

  override protected def coercion = {
    case o: WdlCallOutputsObject => o
  }
}
