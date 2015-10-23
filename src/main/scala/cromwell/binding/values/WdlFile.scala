package cromwell.binding.values

import cromwell.binding.types.{WdlFileType, WdlType}
import cromwell.binding.values.GlobOrNot.GlobOrNot

import scala.util.{Success, Try}

object WdlFile {
  def appendPathsWithSlashSeparators(path1: String, path2: String) = {
    if (path1.endsWith("/") || path2.startsWith("/")) path1 + path2
    else path1 + "/" + path2
  }
}

case class WdlFile(value: String, globbiness: GlobOrNot = GlobOrNot.NotGlob) extends WdlPrimitive {

  val wdlType: WdlType = WdlFileType

  override def add(rhs: WdlValue): Try[WdlValue] = rhs match {
    case r: WdlString => Success(WdlFile(value + r.value))
    case _ => invalid(s"$value + $rhs")
  }

  override def equals(rhs: WdlValue): Try[WdlBoolean] = rhs match {
    case r: WdlFile => Success(WdlBoolean(value.equals(r.value) && globbiness.equals(r.globbiness)))
    case r: WdlString => Success(WdlBoolean(value.toString.equals(r.value.toString) && globbiness == GlobOrNot.NotGlob))
    case _ => invalid(s"$value == $rhs")
  }

  override def toWdlString = "\"" + value.toString + "\""
  override def valueString = value.toString
}

object GlobOrNot extends Enumeration {
  type GlobOrNot = Value
  val YesGlob, NotGlob = Value
}