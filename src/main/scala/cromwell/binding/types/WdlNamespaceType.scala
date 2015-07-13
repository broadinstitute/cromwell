package cromwell.binding.types

case object WdlNamespaceType extends WdlType {
  override def toWdlString: String = "Namespace"
  override protected def coercion = ???
  override def fromRawString(rawString: String) = ???
}
