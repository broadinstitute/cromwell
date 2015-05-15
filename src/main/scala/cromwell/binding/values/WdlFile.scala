package cromwell.binding.values

import java.io.File
import cromwell.binding.types.WdlFileType

import scala.util.{Success, Try}

case class WdlFile(value: File) extends WdlPrimitive {
  val wdlType = WdlFileType
  override def add(rhs: WdlValue): Try[WdlValue] = {
    rhs match {
      case r:WdlString => Success(WdlFile(new File(value + r.value)))
      case r:WdlFile => Success(WdlFile(new File(value + r.value.toString)))
      case _ => invalid(s"$value + $rhs")
    }
  }
  override def equals(rhs: WdlValue): Try[WdlBoolean] = {
    rhs match {
      case r:WdlFile => Success(WdlBoolean(value.equals(r.value)))
      case r:WdlString => Success(WdlBoolean(value.toString.equals(r.value)))
      case _ => invalid(s"$value == $rhs")
    }
  }
  override def asString = value.getAbsolutePath
}
