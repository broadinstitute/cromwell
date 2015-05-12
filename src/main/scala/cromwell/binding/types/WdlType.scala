package cromwell.binding.types

trait WdlType {
  def toWdlString: String
}
