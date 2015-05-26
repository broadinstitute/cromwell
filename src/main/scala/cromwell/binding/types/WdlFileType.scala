package cromwell.binding.types

case object WdlFileType extends WdlType {
  override def toWdlString: String = "File"
}
