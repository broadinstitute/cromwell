package cromwell.binding.expression

import cromwell.binding.WdlExpressionException
import cromwell.binding.types._
import cromwell.binding.values._

import scala.util.{Failure, Success, Try}

trait WdlStandardLibraryFunctions extends WdlFunctions[WdlValue] {
  private def fail(name: String) = Failure(new UnsupportedOperationException(s"$name() not implemented yet"))

  /**
    * Asserts that the parameter list contains a single parameter which will be interpreted
    * as a File and attempts to read the contents of that file and returns back the contents
    * as a String
    */
  private def readContentsFromSingleFileParameter(params: Seq[Try[WdlValue]]): Try[String] = {
    for {
      singleArgument <- extractSingleArgument(params)
      string = fileContentsToString(singleArgument.valueString)
    } yield string
  }

  private def extractObjects(params: Seq[Try[WdlValue]]): Try[Array[WdlObject]] = for {
    contents <- readContentsFromSingleFileParameter(params)
    wdlObjects <- WdlObject.fromTsv(contents)
  } yield wdlObjects

  protected def stdout(params: Seq[Try[WdlValue]]): Try[WdlFile] = fail("stdout")
  protected def stderr(params: Seq[Try[WdlValue]]): Try[WdlFile] = fail("stderr")

  /**
    * Read all lines from the file referenced by the first parameter and return an Array[String]
    */
  protected def read_lines(params: Seq[Try[WdlValue]]): Try[WdlArray] = {
    for {
      contents <- readContentsFromSingleFileParameter(params)
      lines = contents.split("\n")
    } yield WdlArray(WdlArrayType(WdlStringType), lines map WdlString)
  }

  protected def read_map(params: Seq[Try[WdlValue]]): Try[WdlMap] = {
    for {
      contents <- readContentsFromSingleFileParameter(params)
      wdlMap <- WdlMap.fromTsv(contents)
    } yield wdlMap
  }

  protected def read_object(params: Seq[Try[WdlValue]]): Try[WdlObject] = {
    extractObjects(params) map {
      case array if array.length == 1 => array.head
      case _ => throw new IllegalArgumentException("read_object yields an Object and thus can only read 2-rows TSV files. Try using read_objects instead.")
    }
  }

  def read_objects(params: Seq[Try[WdlValue]]): Try[WdlArray] =
    extractObjects(params) map { WdlArray(WdlArrayType(WdlObjectType), _) }

  /**
    * Try to read a string from the file referenced by the specified `WdlValue`.
    */
  protected def read_string(params: Seq[Try[WdlValue]]): Try[WdlString] =
    readContentsFromSingleFileParameter(params).map(s => WdlString(s.trim))

  /**
    * Read a file in TSV format into an Array[Array[String]]
    */
  protected def read_tsv(params: Seq[Try[WdlValue]]): Try[WdlArray] = {
    for {
      contents <- readContentsFromSingleFileParameter(params)
      wdlArray = WdlArray.fromTsv(contents)
    } yield wdlArray
  }

  protected def read_json(params: Seq[Try[WdlValue]]): Try[WdlValue] = fail("read_json")

  /**
   * Try to read an integer from the file referenced by the specified `WdlValue`.
   */
  protected def read_int(params: Seq[Try[WdlValue]]): Try[WdlInteger] =
    read_string(params) map { s => WdlInteger(s.value.trim.toInt) }

  /**
   * Try to read a float from the file referenced by the specified `WdlValue`.
   */
  protected def read_float(params: Seq[Try[WdlValue]]): Try[WdlFloat] =
    read_string(params) map { s => WdlFloat(s.value.trim.toDouble) }

  /**
   * Try to read a boolean from the file referenced by the specified `WdlValue`.
   */
  protected def read_boolean(params: Seq[Try[WdlValue]]): Try[WdlBoolean] =
    read_string(params) map { s => WdlBoolean(java.lang.Boolean.parseBoolean(s.value.trim.toLowerCase)) }
  protected def write_lines(params: Seq[Try[WdlValue]]): Try[WdlFile] = fail("write_lines")
  protected def write_tsv(params: Seq[Try[WdlValue]]): Try[WdlFile] = fail("write_tsv")
  protected def write_map(params: Seq[Try[WdlValue]]): Try[WdlFile] = fail("write_map")
  protected def write_object(params: Seq[Try[WdlValue]]): Try[WdlFile] = fail("write_object")
  protected def write_objects(params: Seq[Try[WdlValue]]): Try[WdlFile] = fail("write_objects")
  protected def write_json(params: Seq[Try[WdlValue]]): Try[WdlFile] = fail("write_json")
  protected def glob(params: Seq[Try[WdlValue]]): Try[WdlArray] = fail("glob")
}

class NoFunctions extends WdlStandardLibraryFunctions

class WdlStandardLibraryFunctionsType extends WdlFunctions[WdlType] {
  protected def stdout(params: Seq[Try[WdlType]]): Try[WdlType] = Success(WdlFileType)
  protected def stderr(params: Seq[Try[WdlType]]): Try[WdlType] = Success(WdlFileType)
  protected def read_lines(params: Seq[Try[WdlType]]): Try[WdlType] = Success(WdlArrayType(WdlStringType))
  protected def read_tsv(params: Seq[Try[WdlType]]): Try[WdlType] = Success(WdlArrayType(WdlArrayType(WdlStringType)))
  protected def read_map(params: Seq[Try[WdlType]]): Try[WdlType] = Success(WdlMapType(WdlStringType, WdlStringType))
  protected def read_object(params: Seq[Try[WdlType]]): Try[WdlType] = Success(WdlObjectType)
  protected def read_objects(params: Seq[Try[WdlType]]): Try[WdlType] = Success(WdlArrayType(WdlObjectType))
  protected def read_json(params: Seq[Try[WdlType]]): Try[WdlType] = Failure(new WdlExpressionException("Return type of read_json() can't be known statically"))
  protected def read_int(params: Seq[Try[WdlType]]): Try[WdlType] = Success(WdlIntegerType)
  protected def read_string(params: Seq[Try[WdlType]]): Try[WdlType] = Success(WdlStringType)
  protected def read_float(params: Seq[Try[WdlType]]): Try[WdlType] = Success(WdlFloatType)
  protected def read_boolean(params: Seq[Try[WdlType]]): Try[WdlType] = Success(WdlBooleanType)
  protected def write_lines(params: Seq[Try[WdlType]]): Try[WdlType] = Success(WdlFileType)
  protected def write_tsv(params: Seq[Try[WdlType]]): Try[WdlType] = Success(WdlFileType)
  protected def write_map(params: Seq[Try[WdlType]]): Try[WdlType] = Success(WdlFileType)
  protected def write_object(params: Seq[Try[WdlType]]): Try[WdlType] = Success(WdlFileType)
  protected def write_objects(params: Seq[Try[WdlType]]): Try[WdlType] = Success(WdlFileType)
  protected def write_json(params: Seq[Try[WdlType]]): Try[WdlType] = Success(WdlFileType)
  protected def glob(params: Seq[Try[WdlType]]): Try[WdlType] = Success(WdlArrayType(WdlFileType))
}
