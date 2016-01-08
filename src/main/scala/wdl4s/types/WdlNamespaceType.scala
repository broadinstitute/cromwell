package wdl4s.types

import wdl4s.WdlNamespace

case object WdlNamespaceType extends WdlType {
  override def toWdlString: String = "Namespace"

  override protected def coercion = {
    case n: WdlNamespace => n
  }

  override def fromWdlString(rawString: String) = ???
}
