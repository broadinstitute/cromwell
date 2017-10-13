package wdl.types

import wdl.WdlNamespace
import wom.types.WdlType

case object WdlNamespaceType extends WdlType {
  override def toWdlString: String = "Namespace"

  override protected def coercion = {
    case n: WdlNamespace => n
  }
}
