package wom.types

import spray.json.JsString
import wom.values.{WomFile, WomString}

import scala.util.{Success, Try}

case object WomFileType extends WomPrimitiveType {
  val toDisplayString: String = "File"

  override protected def coercion = {
    case s: String => WomFile(s)
    case s: JsString => WomFile(s.value)
    case s: WomString => WomFile(s.valueString)
    case f: WomFile => f
  }

  override def add(rhs: WomType): Try[WomType] = rhs match {
    case WomStringType => Success(WomFileType)
    case WomOptionalType(memberType) => add(memberType)
    case _ => invalid(s"$this + $rhs")
  }

  override def equals(rhs: WomType): Try[WomType] = rhs match {
    case WomFileType => Success(WomBooleanType)
    case WomStringType => Success(WomBooleanType)
    case WomOptionalType(memberType) => equals(memberType)
    case _ => invalid(s"$this == $rhs")
  }
}
