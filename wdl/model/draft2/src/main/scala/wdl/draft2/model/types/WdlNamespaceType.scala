package wdl.draft2.model.types

import wdl.draft2.model.WdlNamespace
import wom.types.WomType

case object WdlNamespaceType extends WomType {
  override def stableName: String = "Namespace"

  override protected def coercion = {
    case n: WdlNamespace => n
  }
}
