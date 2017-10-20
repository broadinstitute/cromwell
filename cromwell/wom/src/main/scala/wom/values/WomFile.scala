package wom.values

import wom.types.{WomFileType, WomType}

import scala.util.{Success, Try}

object WomFile {
  def appendPathsWithSlashSeparators(path1: String, path2: String) = {
    if (path1.endsWith("/") || path2.startsWith("/")) path1 + path2
    else path1 + "/" + path2
  }

  def apply(value: String, isGlob: Boolean = false): WomFile = if (isGlob) WomGlobFile(value) else WomSingleFile(value)
}

sealed trait WomFile extends WomPrimitive {
  val value: String
  val womType: WomType = WomFileType

  def isGlob: Boolean = this match {
    case _: WomGlobFile => true
    case _ => false
  }

  override def add(rhs: WomValue): Try[WomValue] = rhs match {
    case r: WomString => Success(WomFile(value + r.value))
    case r: WomOptionalValue => evaluateIfDefined("+", r, add)
    case _ => invalid(s"$value + $rhs")
  }

  override def equals(rhs: WomValue): Try[WomBoolean] = rhs match {
    case r: WomFile => Success(WomBoolean(value.equals(r.value) && isGlob.equals(r.isGlob)))
    case r: WomString => Success(WomBoolean(value.toString.equals(r.value.toString) && !isGlob))
    case r: WomOptionalValue => evaluateIfDefined("==", r, equals)
    case _ => invalid(s"$value == $rhs")
  }

  override def valueString = value.toString
}

case class WomSingleFile(value: String) extends WomFile {
  override def toWomString = "\"" + value.toString + "\""
}

case class WomGlobFile(value: String) extends WomFile {
  override def toWomString = "glob(\"" + value.toString + "\")"
}

