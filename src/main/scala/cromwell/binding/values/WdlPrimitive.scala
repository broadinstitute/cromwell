package cromwell.binding.values

trait WdlPrimitive extends WdlValue {
  def toWdlString: String
  override def toRawString = toWdlString
}
