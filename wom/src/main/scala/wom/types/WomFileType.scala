package wom.types

import spray.json.JsString
import wom.values._

import scala.util.{Success, Try}

sealed trait WomFileType extends WomType {
  override def stableName: String = s"$getClass"

  override protected def coercion: PartialFunction[Any, WomFile] = PartialFunction.empty
}

sealed trait WomPrimitiveFileType extends WomFileType with WomPrimitiveType

case object WomUnlistedDirectoryType extends WomPrimitiveFileType {
  override val stableName: String = "Directory"

  override protected def coercion: PartialFunction[Any, WomUnlistedDirectory] = {
    case s: String => WomUnlistedDirectory(s.trim)
    case s: JsString => WomUnlistedDirectory(s.value.trim)
    case s: WomString => WomUnlistedDirectory(s.valueString.trim)
    case d: WomUnlistedDirectory => d
  }

  override def equalsType(rhs: WomType): Try[WomType] = rhs match {
    case WomUnlistedDirectoryType => Success(WomBooleanType)
    case WomStringType => Success(WomBooleanType)
    case WomOptionalType(memberType) => equalsType(memberType)
    case _ => invalid(s"$this == $rhs")
  }

  override def add(rhs: WomType): Try[WomType] = rhs match {
    case WomStringType => Success(WomUnlistedDirectoryType)
    case WomOptionalType(memberType) => add(memberType)
    case _ => invalid(s"$this + $rhs")
  }
}

case object WomSingleFileType extends WomPrimitiveFileType {
  override val stableName: String = "File"

  override protected def coercion: PartialFunction[Any, WomSingleFile] = {
    case s: String =>
      if (s != "")
        WomSingleFile(s.trim)
      else
        throw new IllegalArgumentException("""Cannot coerce the empty String value "" into a File.""")
    case s: JsString => WomSingleFile(s.value.trim)
    case s: WomString => coercion.apply(s.valueString.trim)
    case f: WomSingleFile => WomSingleFile(f.value.trim)
  }

  override def equalsType(rhs: WomType): Try[WomType] = rhs match {
    case wct:WomCoproductType => wct.typeExists(WomStringType)
    case WomSingleFileType => Success(WomBooleanType)
    case WomStringType => Success(WomBooleanType)
    case WomOptionalType(memberType) => equalsType(memberType)
    case _ => invalid(s"$this == $rhs")
  }

  override def add(rhs: WomType): Try[WomType] = rhs match {
    case WomStringType => Success(WomSingleFileType)
    case WomOptionalType(memberType) => add(memberType)
    case _ => invalid(s"$this + $rhs")
  }
}

case object WomGlobFileType extends WomPrimitiveFileType {
  override val stableName: String = "Glob"

  override protected def coercion: PartialFunction[Any, WomGlobFile] = {
    case s: String => WomGlobFile(s.trim)
    case s: JsString => WomGlobFile(s.value.trim)
    case s: WomString => WomGlobFile(s.valueString.trim)
    case f: WomGlobFile => f
  }

  override def equalsType(rhs: WomType): Try[WomType] = rhs match {
    case wct:WomCoproductType => wct.typeExists(WomStringType)
    case WomGlobFileType => Success(WomBooleanType)
    case WomStringType => Success(WomBooleanType)
    case WomOptionalType(memberType) => equalsType(memberType)
    case _ => invalid(s"$this == $rhs")
  }

  override def add(rhs: WomType): Try[WomType] = rhs match {
    case WomStringType => Success(WomGlobFileType)
    case WomOptionalType(memberType) => add(memberType)
    case _ => invalid(s"$this + $rhs")
  }
}

case object WomMaybeListedDirectoryType extends WomFileType {
  override protected def coercion: PartialFunction[Any, WomMaybeListedDirectory] = {
    case s: String => WomMaybeListedDirectory(s.trim)
    case s: JsString => WomMaybeListedDirectory(s.value.trim)
    case s: WomString => WomMaybeListedDirectory(s.valueString.trim)
    case d: WomUnlistedDirectory => WomMaybeListedDirectory(d.value)
    case d: WomMaybeListedDirectory => d
  }
}

case object WomMaybePopulatedFileType extends WomFileType {
  override protected def coercion: PartialFunction[Any, WomMaybePopulatedFile] = {
    case s: String => WomMaybePopulatedFile(s.trim)
    case s: JsString => WomMaybePopulatedFile(s.value.trim)
    case s: WomString => WomMaybePopulatedFile(s.valueString.trim)
    case f: WomSingleFile => WomMaybePopulatedFile(f.value)
    case f: WomMaybePopulatedFile => f
  }
}
