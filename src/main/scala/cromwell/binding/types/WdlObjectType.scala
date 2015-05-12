package cromwell.binding.types

import cromwell.binding.WdlObject

case object WdlObjectType extends WdlType {
  override def toWdlString: String = "Object"
}
