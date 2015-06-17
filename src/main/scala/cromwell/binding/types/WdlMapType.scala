package cromwell.binding.types

case class WdlMapType(keyType: WdlPrimitiveType, valueType: WdlType) extends WdlType {
  val toWdlString: String = s"Map[${keyType.toWdlString}, ${valueType.toWdlString}]"
  override protected def coercion = ???
  override def fromRawString(s: String) = ???
}
