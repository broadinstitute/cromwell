package cromwell.binding.values

import cromwell.binding.Call
import cromwell.binding.types.{WdlCallOutputsObjectType, WdlObjectType}

trait WdlObjectLike {
  def value: Map[String, WdlValue]
}

case class WdlObject(value: Map[String, WdlValue]) extends WdlValue with WdlObjectLike {
  val wdlType = WdlObjectType
}

case class WdlCallOutputsObject(call: Call, outputs: Map[String, WdlValue]) extends WdlValue with WdlObjectLike {
  val wdlType = WdlCallOutputsObjectType(call)
  val value = outputs
}
