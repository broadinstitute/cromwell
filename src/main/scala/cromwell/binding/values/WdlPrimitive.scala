package cromwell.binding.values

import scala.util.Try

trait WdlPrimitive extends WdlValue {
  def asString: String
  override def toRawString = Try(asString)
}
