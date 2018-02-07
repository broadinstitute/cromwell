package wdl.draft2.types

import wdl.draft2.WdlNamespace
import wom.types.WomType

case object WdlNamespaceType extends WomType {
  override def toDisplayString: String = "Namespace"

  override protected def coercion = {
    case n: WdlNamespace => n
  }
}
