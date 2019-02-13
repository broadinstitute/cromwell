package wdl.draft2.model.types

import wdl.draft2.model.WdlCall
import wdl.draft2.model.values.WdlCallOutputsObject
import wom.types.WomObjectTypeLike

case class WdlCallOutputsObjectType(call: WdlCall) extends WomObjectTypeLike {
  val stableName: String = "Object"

  override protected def coercion = {
    case o: WdlCallOutputsObject => o
  }
}
