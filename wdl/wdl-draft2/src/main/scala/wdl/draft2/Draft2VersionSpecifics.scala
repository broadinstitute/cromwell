package wdl.draft2

import wom.MemorySize
import wdl.versioning.WdlVersionSpecifics
import wdl4s.parser.MemoryUnit
import wom.types.{WomOptionalType, WomSingleFileType, WomType}
import wom.values.{WomFloat, WomOptionalValue, WomSingleFile, WomValue}

import scala.util.{Failure, Success, Try}

case object Draft2VersionSpecifics extends WdlVersionSpecifics {
  override def sizeFunctionOverride(params: Seq[Try[WomValue]], sizeFunc: WomValue => Try[Double]): Try[WomFloat] = {
    // Inner function: get the memory unit from the second (optional) parameter
    def toUnit(womValue: Try[WomValue]) = womValue flatMap { unit => Try(MemoryUnit.fromSuffix(unit.valueString)) }

    // Inner function: is this a file type, or an optional containing a file type?
    def isOptionalOfFileType(womType: WomType): Boolean = womType match {
      case f if WomSingleFileType.isCoerceableFrom(f) => true
      case WomOptionalType(inner) => isOptionalOfFileType(inner)
      case _ => false
    }

    // Inner function: Get the file size, allowing for unpacking of optionals
    def optionalSafeFileSize(value: WomValue): Try[Double] = value match {
      case f if f.isInstanceOf[WomSingleFile] || WomSingleFileType.isCoerceableFrom(f.womType) => sizeFunc(f)
      case WomOptionalValue(_, Some(o)) => optionalSafeFileSize(o)
      case WomOptionalValue(f, None) if isOptionalOfFileType(f) => Success(0d)
      case _ => Failure(new Exception(s"The 'size' method expects a 'File' or 'File?' argument but instead got ${value.womType.toDisplayString}."))
    }

    def fileSize(womValue: Try[WomValue], convertTo: Try[MemoryUnit] = Success(MemoryUnit.Bytes)): Try[Double] = {
      for {
        value <- womValue
        unit <- convertTo
        fileSize <- optionalSafeFileSize(value)
      } yield MemorySize(fileSize, MemoryUnit.Bytes).to(unit).amount
    }

    params match {
      case _ if params.length == 1 => fileSize(params.head) map WomFloat.apply
      case _ if params.length == 2 => fileSize(params.head, toUnit(params.tail.head)) map WomFloat.apply
      case _ => Failure(new UnsupportedOperationException(s"Expected one or two parameters but got ${params.length} instead."))
    }
  }
}
