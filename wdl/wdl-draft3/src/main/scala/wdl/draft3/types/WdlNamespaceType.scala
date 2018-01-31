package wdl.draft3.types

import wdl.draft3.WdlNamespace
import wom.types.WomType

case object WdlNamespaceType extends WomType {
  override def toDisplayString: String = "Namespace"

  override protected def coercion = {
    case n: WdlNamespace => n
  }
}
