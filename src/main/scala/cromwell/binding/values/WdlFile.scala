package cromwell.binding.values

import cromwell.binding.types.{WdlFileType, WdlType}
import cromwell.binding.{FileHasher, SymbolHash}

import scala.util.{Success, Try}

object WdlFile {
  def appendPathsWithSlashSeparators(path1: String, path2: String) = {
    if (path1.endsWith("/") || path2.startsWith("/")) path1 + path2
    else path1 + "/" + path2
  }

  def apply(value: String, isGlob: Boolean = false): WdlFile = if (isGlob) WdlGlobFile(value) else WdlSingleFile(value)
}

sealed trait WdlFile extends WdlPrimitive {
  val value: String

  val wdlType: WdlType = WdlFileType

  def isGlob: Boolean = this match {
    case _: WdlGlobFile => true
    case _ => false
  }

  override def add(rhs: WdlValue): Try[WdlValue] = rhs match {
    case r: WdlString => Success(WdlFile(value + r.value))
    case _ => invalid(s"$value + $rhs")
  }

  override def equals(rhs: WdlValue): Try[WdlBoolean] = rhs match {
    case r: WdlFile => Success(WdlBoolean(value.equals(r.value) && isGlob.equals(r.isGlob)))
    case r: WdlString => Success(WdlBoolean(value.toString.equals(r.value.toString) && !isGlob))
    case _ => invalid(s"$value == $rhs")
  }
  override def valueString = value.toString

  override def getHash(implicit hasher: FileHasher): SymbolHash = hasher(this)
}

case class WdlSingleFile(value: String) extends WdlFile {
  override def toWdlString = "\"" + value.toString + "\""
}

case class WdlGlobFile(value: String) extends WdlFile {
  override def toWdlString = "glob(\"" + value.toString + "\")"
}

