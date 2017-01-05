package cromwell.backend.wdl

import cromwell.backend.MemorySize
import cromwell.core.path.PathFactory
import wdl4s.expression.WdlStandardLibraryFunctions
import wdl4s.parser.MemoryUnit
import wdl4s.types.{WdlArrayType, WdlFileType, WdlObjectType, WdlStringType}
import wdl4s.values._

import scala.util.{Failure, Success, Try}

trait ReadLikeFunctions extends PathFactory { this: WdlStandardLibraryFunctions =>
  import better.files._

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

  private def extractObjects(functionName: String, params: Seq[Try[WdlValue]]): Try[Array[WdlObject]] = for {
    contents <- readContentsFromSingleFileParameter(functionName, params)
    wdlObjects <- WdlObject.fromTsv(contents)
  } yield wdlObjects

  override def readFile(path: String): String = File(buildPath(path)).contentAsString

  /**
    * Read all lines from the file referenced by the first parameter and return an Array[String]
    */
  override def read_lines(params: Seq[Try[WdlValue]]): Try[WdlArray] = {
    for {
      contents <- readContentsFromSingleFileParameter("read_lines", params)
      lines = contents.split("\n")
    } yield WdlArray(WdlArrayType(WdlStringType), lines map WdlString)
  }

  override def read_map(params: Seq[Try[WdlValue]]): Try[WdlMap] = {
    for {
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
  override def read_string(params: Seq[Try[WdlValue]]): Try[WdlString] = readContentsFromSingleFileParameter("read_string", params).map(s => WdlString(s.trim))

  /**
    * Read a file in TSV format into an Array[Array[String]]
    */
  override def read_tsv(params: Seq[Try[WdlValue]]): Try[WdlArray] = {
    for {
      contents <- readContentsFromSingleFileParameter("read_tsv", params)
      wdlArray = WdlArray.fromTsv(contents)
    } yield wdlArray
  }

  /**
    * Try to read an integer from the file referenced by the specified `WdlValue`.
    */
  override def read_int(params: Seq[Try[WdlValue]]): Try[WdlInteger] = read_string(params) map { s => WdlInteger(s.value.trim.toInt) }

  /**
    * Try to read a float from the file referenced by the specified `WdlValue`.
    */
  override def read_float(params: Seq[Try[WdlValue]]): Try[WdlFloat] = read_string(params) map { s => WdlFloat(s.value.trim.toDouble) }

  /**
    * Try to read a boolean from the file referenced by the specified `WdlValue`.
    */
  override def read_boolean(params: Seq[Try[WdlValue]]): Try[WdlBoolean] =
    read_string(params) map { s => WdlBoolean(java.lang.Boolean.parseBoolean(s.value.trim.toLowerCase)) }

  override def size(params: Seq[Try[WdlValue]]): Try[WdlFloat] = {
    def toUnit(wdlValue: Try[WdlValue]) = wdlValue flatMap { unit => Try(MemoryUnit.fromSuffix(unit.valueString)) }

    def fileSize(wdlValue: Try[WdlValue], convertTo: Try[MemoryUnit] = Success(MemoryUnit.Bytes)) = {
      for {
        value <- wdlValue
        unit <- convertTo
      } yield MemorySize(File(buildPath(value.valueString)).size.toDouble, MemoryUnit.Bytes).to(unit).amount
    }

    params match {
      case oneArg if params.length == 1 => fileSize(params.head) map WdlFloat.apply
      case twoArgs if params.length == 2 => fileSize(params.head, toUnit(params(1))) map WdlFloat.apply
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
