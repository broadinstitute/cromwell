package cromwell.binding.types

case object WdlBooleanType extends WdlType {
  override def toWdlString: String = "Boolean"
}

