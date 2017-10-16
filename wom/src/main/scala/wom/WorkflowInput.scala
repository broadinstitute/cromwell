package wom

import wom.core.FullyQualifiedName
import wom.types.{WomOptionalType, WomType}

case class WorkflowInput(fqn: FullyQualifiedName, womType: WomType) {
  val optional = womType.isInstanceOf[WomOptionalType]
}
