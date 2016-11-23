package wdl4s

import wdl4s.types.{WdlOptionalType, WdlType}

case class WorkflowInput(fqn: FullyQualifiedName, wdlType: WdlType) {
  val optional = wdlType.isInstanceOf[WdlOptionalType]
}
