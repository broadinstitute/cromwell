package cromwell.binding.types

case object WdlObjectType extends WdlType {
  override def toWdlString: String = "Object"
}
