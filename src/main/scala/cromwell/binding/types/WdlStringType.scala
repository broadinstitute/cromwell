package cromwell.binding.types

case object WdlStringType extends WdlType {
  override def toWdlString: String = "String"
}
