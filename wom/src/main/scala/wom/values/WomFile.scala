package wom.values

import wom.types.{WomFileType, WomGlobFileType, WomSingleFileType, WomType}

import scala.util.{Success, Try}

sealed trait WomFile extends WomPrimitive {
  val value: String

  override def valueString = value

  // TODO: WOM: WOMFILE: When WDL supports directories this check may need more work (refactoring to subclasses?)
  override def equals(rhs: WomValue): Try[WomBoolean] = rhs match {
    case r: WomFile => Success(WomBoolean(value.equals(r.value) && womType.equals(r.womType)))
    case r: WomString => Success(WomBoolean(value.toString.equals(r.value.toString) && womType == WomSingleFileType))
    case r: WomOptionalValue => evaluateIfDefined("==", r, equals)
    case _ => invalid(s"$value == $rhs")
  }
}

object WomFile {
  def apply(fileType: WomFileType, value: String) = {
    fileType match {
      case WomSingleFileType => WomSingleFile(value)
      case WomGlobFileType => WomGlobFile(value)
    }
  }

  def mapFile(womFile: WomFile, f: String => String): WomFile = {
    womFile match {
      case womSingleFile: WomSingleFile => mapSingleFile(womSingleFile, f)
      case womGlobFile: WomGlobFile => mapGlobFile(womGlobFile, f)
    }
  }

  def mapSingleFile(womFile: WomSingleFile, f: String => String): WomSingleFile = {
    WomSingleFile(f(womFile.value))
  }

  def mapGlobFile(womFile: WomGlobFile, f: String => String): WomGlobFile = {
    WomGlobFile(f(womFile.value))
  }
}

case class WomSingleFile(value: String) extends WomFile {
  override val womType: WomType = WomSingleFileType

  override def toWomString = "\"" + value.toString + "\""

  override def add(rhs: WomValue): Try[WomValue] = rhs match {
    case r: WomString => Success(WomSingleFile(value + r.value))
    case r: WomOptionalValue => evaluateIfDefined("+", r, add)
    case _ => invalid(s"$value + $rhs")
  }
}

case class WomGlobFile(value: String) extends WomFile {
  override val womType: WomType = WomGlobFileType

  override def toWomString = "glob(\"" + value.toString + "\")"

  override def add(rhs: WomValue): Try[WomValue] = rhs match {
    case r: WomString => Success(WomGlobFile(value + r.value))
    case r: WomOptionalValue => evaluateIfDefined("+", r, add)
    case _ => invalid(s"$value + $rhs")
  }
}
