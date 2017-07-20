package cromwell.backend.wdl

import cromwell.backend.MemorySize
import cromwell.core.path.PathFactory
import wdl4s.wdl.expression.WdlStandardLibraryFunctions
import wdl4s.parser.MemoryUnit
import wdl4s.wdl.types._
import wdl4s.wdl.values._

import scala.util.{Failure, Success, Try}

trait ReadLikeFunctions extends PathFactory { this: WdlStandardLibraryFunctions =>

  val fileSizeLimitationConfig =  FileSizeLimitationConfig.fileSizeLimitationConfig
  import fileSizeLimitationConfig._

  def fileSize: WdlValue=> Try[Long] = 
    w => Try(buildPath(w.valueString).size)

  /**
    * Asserts that the parameter list contains a single parameter which will be interpreted
    * as a File and attempts to read the contents of that file
    */
  private def readContentsFromSingleFileParameter(functionName: String, params: Seq[Try[WdlValue]]): Try[String] = {
    for {
      singleArgument <- extractSingleArgument(functionName, params)
      string = fileContentsToString(singleArgument.valueString)
    } yield string
  }

  private def extractObjects(functionName: String, params: Seq[Try[WdlValue]]): Try[Array[WdlObject]] =
    for {
      _ <- validateFileSizeIsWithinLimits("read_object", params, readObjectLimit)
      contents <- readContentsFromSingleFileParameter(functionName, params)
      wdlObjects <- WdlObject.fromTsv(contents)
    } yield wdlObjects

  override def readFile(path: String): String = buildPath(path).contentAsString

  /**
    * Read all lines from the file referenced by the first parameter and return an Array[String]
    */
  override def read_lines(params: Seq[Try[WdlValue]]): Try[WdlArray] = {
    for {
      _ <- validateFileSizeIsWithinLimits("read_lines", params, readLinesLimit)
      contents <- readContentsFromSingleFileParameter("read_lines", params)
      lines = contents.split("\n")
    } yield WdlArray(WdlArrayType(WdlStringType), lines map WdlString)
  }

  def validateFileSizeIsWithinLimits(functionName: String, params: Seq[Try[WdlValue]], limit: Int): Try[Unit] = 
    for {
      fileName <- extractSingleArgument(functionName, params)
      fileSize <- fileSize(fileName)
      _ = if (fileSize > limit) {
        val errorMsg = s"Use of $fileName failed because the file was too big ($fileSize bytes when only files of up to $limit bytes are permissible"
        throw FileSizeTooBig(errorMsg)
      }
    } yield ()

  override def read_map(params: Seq[Try[WdlValue]]): Try[WdlMap] = {
    for {
      _ <- validateFileSizeIsWithinLimits("read_map", params, readMapLimit)
      contents <- readContentsFromSingleFileParameter("read_map", params)
      wdlMap <- WdlMap.fromTsv(contents)
    } yield wdlMap
  }

  override def read_object(params: Seq[Try[WdlValue]]): Try[WdlObject] = {
    extractObjects("read_object", params) map {
      case array if array.length == 1 => array.head
      case _ => throw new IllegalArgumentException("read_object yields an Object and thus can only read 2-rows TSV files. Try using read_objects instead.")
    }
  }

  override def read_objects(params: Seq[Try[WdlValue]]): Try[WdlArray] = extractObjects("read_objects", params) map { WdlArray(WdlArrayType(WdlObjectType), _) }

  /**
    * Try to read a string from the file referenced by the specified `WdlValue`.
    */
  override def read_string(params: Seq[Try[WdlValue]]): Try[WdlString] =
    for {
      _ <- validateFileSizeIsWithinLimits("read_string", params, readStringLimit)
      string <- readContentsFromSingleFileParameter("read_string", params)
    } yield WdlString(string.trim)

  /**
    * Read a file in TSV format into an Array[Array[String]]
    */
  override def read_tsv(params: Seq[Try[WdlValue]]): Try[WdlArray] = {
    for {
      _ <- validateFileSizeIsWithinLimits("read_tsv", params, readTsvLimit)
      contents <- readContentsFromSingleFileParameter("read_tsv", params)
      wdlArray = WdlArray.fromTsv(contents)
    } yield wdlArray
  }

  /**
    * Try to read an integer from the file referenced by the specified `WdlValue`.
    */
  override def read_int(params: Seq[Try[WdlValue]]): Try[WdlInteger] = 
    for {
      _ <- validateFileSizeIsWithinLimits("read_int", params, readIntLimit)
      r <- read_string(params) map { s => WdlInteger(s.value.trim.toInt) }
    } yield r

  /**
    * Try to read a float from the file referenced by the specified `WdlValue`.
    */
  override def read_float(params: Seq[Try[WdlValue]]): Try[WdlFloat] =
    for {
      _ <- validateFileSizeIsWithinLimits("read_float", params, readFloatLimit)
      s <- read_string(params)
    } yield WdlFloat(s.value.trim.toDouble)

  /**
    * Try to read a boolean from the file referenced by the specified `WdlValue`.
    */
  override def read_boolean(params: Seq[Try[WdlValue]]): Try[WdlBoolean] =
    read_string(params) map { s => WdlBoolean(java.lang.Boolean.parseBoolean(s.value.trim.toLowerCase)) }

  protected def size(file: WdlValue): Try[Double] = Try(buildPath(file.valueString).size.toDouble)

  /**
    * Gets the size of a file.
    *
    * @param params First parameter must be a File or File? or coerceable to one. The second is an optional string containing the size unit (eg "MB", "GiB")
    */
  override def size(params: Seq[Try[WdlValue]]): Try[WdlFloat] = {
    // Inner function: get the memory unit from the second (optional) parameter
    def toUnit(wdlValue: Try[WdlValue]) = wdlValue flatMap { unit => Try(MemoryUnit.fromSuffix(unit.valueString)) }

    // Inner function: is this a file type, or an optional containing a file type?
    def isOptionalOfFileType(wdlType: WdlType): Boolean = wdlType match {
      case f if WdlFileType.isCoerceableFrom(f) => true
      case WdlOptionalType(inner) => isOptionalOfFileType(inner)
      case _ => false
    }

    // Inner function: Get the file size, allowing for unpacking of optionals
    def optionalSafeFileSize(value: WdlValue): Try[Double] = value match {
      case f if f.isInstanceOf[WdlFile] || WdlFileType.isCoerceableFrom(f.wdlType) => size(f)
      case WdlOptionalValue(_, Some(o)) => optionalSafeFileSize(o)
      case WdlOptionalValue(f, None) if isOptionalOfFileType(f) => Success(0d)
      case _ => Failure(new Exception(s"The 'size' method expects a 'File' or 'File?' argument but instead got ${value.wdlType.toWdlString}."))
    }

    // Inner function: get the file size and convert into the requested memory unit
    def fileSize(wdlValue: Try[WdlValue], convertTo: Try[MemoryUnit] = Success(MemoryUnit.Bytes)) = {
      for {
        value <- wdlValue
        unit <- convertTo
        fileSize <- optionalSafeFileSize(value)
      } yield MemorySize(fileSize, MemoryUnit.Bytes).to(unit).amount
    }

    params match {
      case _ if params.length == 1 => fileSize(params.head) map WdlFloat.apply
      case _ if params.length == 2 => fileSize(params.head, toUnit(params.tail.head)) map WdlFloat.apply
      case _ => Failure(new UnsupportedOperationException(s"Expected one or two parameters but got ${params.length} instead."))
    }
  }

  override def glob(params: Seq[Try[WdlValue]]): Try[WdlArray] = {
    for {
      singleArgument <- extractSingleArgument("glob", params)
      globVal = singleArgument.valueString
      files = glob(globPath(globVal), globVal)
      wdlFiles = files map { WdlFile(_, isGlob = false) }
    } yield WdlArray(WdlArrayType(WdlFileType), wdlFiles)
  }

  override def read_json(params: Seq[Try[WdlValue]]): Try[WdlValue] = Failure(new NotImplementedError(s"read_json() not implemented yet"))

}
