package cromwell.binding.values

import cromwell.binding.types.WdlObjectType

case class WdlObject(value: Map[String, WdlValue]) extends WdlValue {
  val wdlType = WdlObjectType

  override def toRawString = ???
  override def toWdlString = ???
}
