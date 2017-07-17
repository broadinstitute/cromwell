package wdl4s.wdl

import wdl4s.wdl.types.{WdlOptionalType, WdlType}
import wdl4s.wom.callable.Callable

case class WorkflowInput(fqn: FullyQualifiedName, wdlType: WdlType) {
  val optional = wdlType.isInstanceOf[WdlOptionalType]
  def toWom = Callable.RequiredInputDefinition(fqn, wdlType)
}
