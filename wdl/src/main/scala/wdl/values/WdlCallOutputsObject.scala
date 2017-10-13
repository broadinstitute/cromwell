package wdl.values

import wdl.WdlCall
import wdl.types.WdlCallOutputsObjectType
import wom.values.{WdlObjectLike, WdlValue}

case class WdlCallOutputsObject(call: WdlCall, outputs: Map[String, WdlValue]) extends WdlValue with WdlObjectLike {
  val wdlType = WdlCallOutputsObjectType(call)
  val value = outputs

  override def copyWith(values: Map[String, WdlValue]) = copy(call, values)
}
