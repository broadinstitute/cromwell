package cromwell.binding.types

case object WdlIntegerType extends WdlType {
  override def toWdlString: String = "Int"
}