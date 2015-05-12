package cromwell.binding.types

case object WdlFloatType extends WdlType {
  override def toWdlString: String = "Float"
}

