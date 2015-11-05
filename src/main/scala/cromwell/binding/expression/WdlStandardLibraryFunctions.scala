package cromwell.binding.expression

import cromwell.binding.types._
import cromwell.binding.values._
import cromwell.binding.{IoInterface, TsvSerializable, WdlExpressionException}

import scala.language.postfixOps
import scala.util.{Failure, Success, Try}

trait WdlStandardLibraryFunctions extends WdlFunctions[WdlValue] {
  def interface: IoInterface

  private def fail(name: String) = Failure(new NotImplementedError(s"$name() not implemented yet"))

  protected def fileContentsToString(path: String): String = interface.readFile(path)

  private def writeContent(baseName: String, content: String): Try[WdlFile] = {
    Try(WdlFile(interface.writeTempFile(tempFilePath, s"$baseName.", ".tmp", content)))
  }

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

  private def writeToTsv(params: Seq[Try[WdlValue]], wdlClass: Class[_ <: WdlValue with TsvSerializable]) = {
    for {
      singleArgument <- extractSingleArgument(params)
      downcast <- Try(wdlClass.cast(singleArgument))
      tsvSerialized <- downcast.tsvSerialize
      file <- writeContent(wdlClass.getSimpleName.toLowerCase, tsvSerialized)
    } yield file
  }

  protected def write_lines(params: Seq[Try[WdlValue]]): Try[WdlFile] = writeToTsv(params, classOf[WdlArray])

  protected def write_map(params: Seq[Try[WdlValue]]): Try[WdlFile] = writeToTsv(params, classOf[WdlMap])

  protected def write_object(params: Seq[Try[WdlValue]]): Try[WdlFile] = writeToTsv(params, classOf[WdlObject])

  protected def write_objects(params: Seq[Try[WdlValue]]): Try[WdlFile] = writeToTsv(params, classOf[WdlArray])

  protected def glob(params: Seq[Try[WdlValue]]): Try[WdlArray] = {
    for {
      singleArgument <- extractSingleArgument(params)
      globVal = singleArgument.valueString
      files = interface.glob(globPath(globVal), globVal)
      wdlFiles = files map { WdlFile(_, isGlob = false) }
    } yield WdlArray(WdlArrayType(WdlFileType), wdlFiles toSeq)
  }

  protected def write_tsv(params: Seq[Try[WdlValue]]): Try[WdlFile] = fail("write_tsv")
  protected def write_json(params: Seq[Try[WdlValue]]): Try[WdlFile] = fail("write_json")
}

class NoFunctions extends WdlStandardLibraryFunctions with IoInterface {
  override def readFile(path: String): String = throw new NotImplementedError()
  override def writeFile(path: String, content: String): Unit = throw new NotImplementedError()
  override def listContents(path: String): Iterable[String] = throw new NotImplementedError()
  override def glob(path: String, pattern: String): Seq[String] = throw new NotImplementedError()
  override def writeTempFile(path: String, prefix: String, suffix: String, content: String): String = throw new NotImplementedError()
  override def exists(path: String): Boolean = throw new NotImplementedError()
  override def interface: IoInterface = throw new NotImplementedError()
  override def isValidPath(path: String): Boolean = throw new NotImplementedError()
  override def copy(from: String, to: String): Unit = throw new NotImplementedError()
  override def hash(path: String): String = throw new NotImplementedError()
}

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
