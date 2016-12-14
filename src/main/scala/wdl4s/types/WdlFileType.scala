package wdl4s.types

import wdl4s.values.{WdlFile, WdlString}
import spray.json.JsString

import scala.util.{Success, Try}

case object WdlFileType extends WdlPrimitiveType {
  val toWdlString: String = "File"

  override protected def coercion = {
    case s: String => WdlFile(s)
    case s: JsString => WdlFile(s.value)
    case s: WdlString => WdlFile(s.valueString)
    case f: WdlFile => f
  }

  override def add(rhs: WdlType): Try[WdlType] = rhs match {
    case WdlStringType => Success(WdlFileType)
    case WdlOptionalType(memberType) => add(memberType)
    case _ => invalid(s"$this + $rhs")
  }

  override def equals(rhs: WdlType): Try[WdlType] = rhs match {
    case WdlFileType => Success(WdlBooleanType)
    case WdlStringType => Success(WdlBooleanType)
    case WdlOptionalType(memberType) => equals(memberType)
    case _ => invalid(s"$this == $rhs")
  }
}
