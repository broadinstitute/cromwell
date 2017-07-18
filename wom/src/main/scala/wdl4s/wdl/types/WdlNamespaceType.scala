package wdl4s.wdl.types

import wdl4s.wdl.WdlNamespace

case object WdlNamespaceType extends WdlType {
  override def toWdlString: String = "Namespace"

  override protected def coercion = {
    case n: WdlNamespace => n
  }

  override def fromWdlString(rawString: String) = ???
}
