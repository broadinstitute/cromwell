package cromwell.binding.values

trait WdlPrimitive extends WdlValue {
  def asString: String
}
