package cromwell.binding.types

import cromwell.binding.Call
import cromwell.binding.values.{WdlCallOutputsObject, WdlObject}

case object WdlObjectType extends WdlType {
  val toWdlString: String = "Object"

  override protected def coercion = {
    case o: WdlObject => o
  }

  override def fromWdlString(rawString: String) = ???
}

case class WdlCallOutputsObjectType(call: Call) extends WdlType {
  val toWdlString: String = "Object"

  override protected def coercion = {
    case o: WdlCallOutputsObject => o
  }
}
