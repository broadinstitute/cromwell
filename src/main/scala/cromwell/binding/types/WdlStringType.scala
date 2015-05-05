package cromwell.binding.types

case object WdlStringType extends WdlType {
    def isCompatible(value: Any) = value.isInstanceOf[String]
    override def toString: String = "string"
}
