package wdl.draft2.model.types

import wdl.draft2.model.WdlNamespace
import wom.types.WomType
import wom.values.WomValue

case object WdlNamespaceType extends WomType {
  override def toDisplayString: String = "Namespace"

  override def coercion: PartialFunction[Any, WomValue] = {
    case n: WdlNamespace => n
  }
}
