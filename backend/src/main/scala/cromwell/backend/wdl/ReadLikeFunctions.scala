package cromwell.backend.wdl

import cromwell.backend.MemorySize
import cromwell.core.path.PathFactory
import wdl4s.parser.MemoryUnit
import wom.expression.IoFunctionSet
import wom.types._
import wom.values._

import scala.concurrent.Future
import scala.util.{Failure, Success, Try}

trait ReadLikeFunctions extends PathFactory with IoFunctionSet {

  // TODO WOM: https://github.com/broadinstitute/cromwell/issues/2611
  val fileSizeLimitationConfig =  FileSizeLimitationConfig.fileSizeLimitationConfig

  def fileSize: WomValue=> Try[Long] =
    w => Try(buildPath(w.valueString).size)

  // TODO WOM: https://github.com/broadinstitute/cromwell/issues/2612
  override def readFile(path: String): Future[String] = Future.successful(buildPath(path).contentAsString)

  protected def size(file: WomValue): Try[Double] = Try(buildPath(file.valueString).size.toDouble)

  /**
    * Gets the size of a file.
    *
    * @param params First parameter must be a File or File? or coerceable to one. The second is an optional string containing the size unit (eg "MB", "GiB")
    */
  override def size(params: Seq[Try[WomValue]]): Try[WomFloat] = {
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
      case f if f.isInstanceOf[WomSingleFile] || WomSingleFileType.isCoerceableFrom(f.womType) => size(f)
      case WomOptionalValue(_, Some(o)) => optionalSafeFileSize(o)
      case WomOptionalValue(f, None) if isOptionalOfFileType(f) => Success(0d)
      case _ => Failure(new Exception(s"The 'size' method expects a 'File' or 'File?' argument but instead got ${value.womType.toDisplayString}."))
    }

    // Inner function: get the file size and convert into the requested memory unit
    def fileSize(womValue: Try[WomValue], convertTo: Try[MemoryUnit] = Success(MemoryUnit.Bytes)) = {
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
