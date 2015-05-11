package cromwell.binding.types

case object WdlFloatType extends WdlType {
  def isCompatible(value: Any) = value.isInstanceOf[Float]

  override def toString: String = "Float"
}

