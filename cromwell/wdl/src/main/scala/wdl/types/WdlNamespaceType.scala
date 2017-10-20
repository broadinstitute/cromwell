package wdl.types

import wdl.WdlNamespace
import wom.types.WomType

case object WdlNamespaceType extends WomType {
  override def toDisplayString: String = "Namespace"

  override protected def coercion = {
    case n: WdlNamespace => n
  }
}
