package cromwell.binding.types

import java.io.File

case object WdlFileType extends WdlType {
  def isCompatible(value: Any) = value.isInstanceOf[File]

  override def toString: String = "File"
}
