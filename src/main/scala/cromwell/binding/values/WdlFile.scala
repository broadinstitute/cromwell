package cromwell.binding.values

import cromwell.binding.types.{WdlType, WdlFileType}

import scala.util.{Success, Try}

object WdlFile {
  def appendPathsWithSlashSeparators(path1: String, path2: String) = {
    if(path1.endsWith("/") || path2.startsWith("/")) path1 + path2
    else path1 + "/" + path2
  }
}

case class WdlFile(value: String) extends WdlPrimitive {

  val wdlType: WdlType = WdlFileType

  override def add(rhs: WdlValue): Try[WdlValue] = {
    rhs match {
      case r: WdlString => Success(WdlFile( WdlFile.appendPathsWithSlashSeparators(value, r.value)))
      case r: WdlFile => Success(WdlFile(WdlFile.appendPathsWithSlashSeparators(value, r.value)))
      case _ => invalid(s"$value + $rhs")
    }
  }
  override def equals(rhs: WdlValue): Try[WdlBoolean] = {
    rhs match {
      case r: WdlFile => Success(WdlBoolean(value.equals(r.value)))
      case r: WdlString => Success(WdlBoolean(value.toString.equals(r.value)))
      case _ => invalid(s"$value == $rhs")
    }
  }

  override def asString = value
}
