package wdl

import wdl.types.{WdlOptionalType, WdlType}
import wom.callable.Callable

case class WorkflowInput(fqn: FullyQualifiedName, wdlType: WdlType) {
  val optional = wdlType.isInstanceOf[WdlOptionalType]
  def toWom = Callable.RequiredInputDefinition(fqn, wdlType)
}
