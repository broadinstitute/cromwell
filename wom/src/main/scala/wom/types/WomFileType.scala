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
    case s: String => WomUnlistedDirectory(s)
    case s: JsString => WomUnlistedDirectory(s.value)
    case s: WomString => WomUnlistedDirectory(s.valueString)
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
        WomSingleFile(s)
      else
        throw new IllegalArgumentException("""Cannot coerce the empty String value "" into a File.""")
    case s: JsString => WomSingleFile(s.value)
    case s: WomString => coercion.apply(s.valueString)
    case f: WomSingleFile => f
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
    case s: String => WomGlobFile(s)
    case s: JsString => WomGlobFile(s.value)
    case s: WomString => WomGlobFile(s.valueString)
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
    case s: String => WomMaybeListedDirectory(s)
    case s: JsString => WomMaybeListedDirectory(s.value)
    case s: WomString => WomMaybeListedDirectory(s.valueString)
    case d: WomUnlistedDirectory => WomMaybeListedDirectory(d.value)
    case d: WomMaybeListedDirectory => d
  }
}

case object WomMaybePopulatedFileType extends WomFileType {
  override protected def coercion: PartialFunction[Any, WomMaybePopulatedFile] = {
    case s: String => WomMaybePopulatedFile(s)
    case s: JsString => WomMaybePopulatedFile(s.value)
    case s: WomString => WomMaybePopulatedFile(s.valueString)
    case f: WomSingleFile => WomMaybePopulatedFile(f.value)
    case f: WomMaybePopulatedFile => f
  }
}
