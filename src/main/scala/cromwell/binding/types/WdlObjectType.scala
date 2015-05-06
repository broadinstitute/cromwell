package cromwell.binding.types

import cromwell.binding.WdlObject

case object WdlObjectType extends WdlType {
  def isCompatible(value: Any) = value.isInstanceOf[WdlObject]

  override def toString: String = "Object"
}
