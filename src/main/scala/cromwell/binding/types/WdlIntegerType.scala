package cromwell.binding.types

case object WdlIntegerType extends WdlType {
  def isCompatible(value: Any) = value.isInstanceOf[Integer]

  override def toString: String = "int"
}