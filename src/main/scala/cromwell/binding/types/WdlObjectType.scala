package cromwell.binding.types

case object WdlObjectType extends WdlType {
  val toWdlString: String = "Object"

  override protected def coercion = ???

  override def fromRawString(rawString: String) = ???
}
