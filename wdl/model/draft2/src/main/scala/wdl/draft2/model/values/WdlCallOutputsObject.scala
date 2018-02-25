package wdl.draft2.model.values

import wdl.draft2.model.WdlCall
import wdl.draft2.model.types.WdlCallOutputsObjectType
import wom.values.{WomObjectLike, WomValue}

case class WdlCallOutputsObject(call: WdlCall, outputs: Map[String, WomValue]) extends WomValue with WomObjectLike {
  val womType = WdlCallOutputsObjectType(call)
  val values = outputs
  lazy val womObjectTypeLike = womType
  override def copyWith(values: Map[String, WomValue]) = copy(call, values)
}
