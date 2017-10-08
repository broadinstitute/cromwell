package wdl

import wdl.types.{WdlOptionalType, WdlType}

case class WorkflowInput(fqn: FullyQualifiedName, wdlType: WdlType) {
  val optional = wdlType.isInstanceOf[WdlOptionalType]
}
