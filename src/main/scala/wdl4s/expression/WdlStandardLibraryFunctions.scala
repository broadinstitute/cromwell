package wdl4s.expression

import wdl4s.types._
import wdl4s.values._
import wdl4s.{TsvSerializable, WdlExpressionException}

import scala.language.postfixOps
import scala.util.{Failure, Success, Try}

trait WdlStandardLibraryFunctions extends WdlFunctions[WdlValue] {
  def fileContentsToString(path: String): String = readFile(path)
  def readFile(path: String): String
  def writeTempFile(path: String, prefix: String, suffix: String, content: String): String
  def stdout(params: Seq[Try[WdlValue]]): Try[WdlFile]
  def stderr(params: Seq[Try[WdlValue]]): Try[WdlFile]
  def glob(path: String, pattern: String): Seq[String]
  def write_tsv(params: Seq[Try[WdlValue]]): Try[WdlFile]
  def write_json(params: Seq[Try[WdlValue]]): Try[WdlFile]
  def size(params: Seq[Try[WdlValue]]): Try[WdlFloat]
  def sub(params: Seq[Try[WdlValue]]): Try[WdlString]

  def read_objects(params: Seq[Try[WdlValue]]): Try[WdlArray] = extractObjects(params) map { WdlArray(WdlArrayType(WdlObjectType), _) }
  def read_string(params: Seq[Try[WdlValue]]): Try[WdlString] = readContentsFromSingleFileParameter(params).map(s => WdlString(s.trim))
  def read_json(params: Seq[Try[WdlValue]]): Try[WdlValue]
  def read_int(params: Seq[Try[WdlValue]]): Try[WdlInteger] = read_string(params) map { s => WdlInteger(s.value.trim.toInt) }
  def read_float(params: Seq[Try[WdlValue]]): Try[WdlFloat] = read_string(params) map { s => WdlFloat(s.value.trim.toDouble) }

  def write_lines(params: Seq[Try[WdlValue]]): Try[WdlFile] = writeToTsv(params, classOf[WdlArray])
  def write_map(params: Seq[Try[WdlValue]]): Try[WdlFile] = writeToTsv(params, classOf[WdlMap])
  def write_object(params: Seq[Try[WdlValue]]): Try[WdlFile] = writeToTsv(params, classOf[WdlObject])
  def write_objects(params: Seq[Try[WdlValue]]): Try[WdlFile] = writeToTsv(params, classOf[WdlArray])

  def read_lines(params: Seq[Try[WdlValue]]): Try[WdlArray] = {
    for {
      contents <- readContentsFromSingleFileParameter(params)
      lines = contents.split("\n")
    } yield WdlArray(WdlArrayType(WdlStringType), lines map WdlString)
  }

  def read_map(params: Seq[Try[WdlValue]]): Try[WdlMap] = {
    for {
      contents <- readContentsFromSingleFileParameter(params)
      wdlMap <- WdlMap.fromTsv(contents)
    } yield wdlMap
  }

  def read_object(params: Seq[Try[WdlValue]]): Try[WdlObject] = {
    extractObjects(params) map {
      case array if array.length == 1 => array.head
      case _ => throw new IllegalArgumentException("read_object yields an Object and thus can only read 2-rows TSV files. Try using read_objects instead.")
    }
  }

  def read_tsv(params: Seq[Try[WdlValue]]): Try[WdlArray] = {
    for {
      contents <- readContentsFromSingleFileParameter(params)
      wdlArray = WdlArray.fromTsv(contents)
    } yield wdlArray
  }

  def read_boolean(params: Seq[Try[WdlValue]]): Try[WdlBoolean] = {
    read_string(params) map { s => WdlBoolean(java.lang.Boolean.parseBoolean(s.value.trim.toLowerCase)) }
  }

  def glob(params: Seq[Try[WdlValue]]): Try[WdlArray] = {
    for {
      singleArgument <- extractSingleArgument(params)
      globVal = singleArgument.valueString
      files = glob(globPath(globVal), globVal)
      wdlFiles = files map { WdlFile(_, isGlob = false) }
    } yield WdlArray(WdlArrayType(WdlFileType), wdlFiles toSeq)
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

  private def writeContent(baseName: String, content: String): Try[WdlFile] = {
    Try(WdlFile(writeTempFile(tempFilePath, s"$baseName.", ".tmp", content)))
  }

  private def writeToTsv(params: Seq[Try[WdlValue]], wdlClass: Class[_ <: WdlValue with TsvSerializable]) = {
    for {
      singleArgument <- extractSingleArgument(params)
      downcast <- Try(wdlClass.cast(singleArgument))
      tsvSerialized <- downcast.tsvSerialize
      file <- writeContent(wdlClass.getSimpleName.toLowerCase, tsvSerialized)
    } yield file
  }
}

class WdlStandardLibraryFunctionsType extends WdlFunctions[WdlType] {
  def stdout(params: Seq[Try[WdlType]]): Try[WdlType] = Success(WdlFileType)
  def stderr(params: Seq[Try[WdlType]]): Try[WdlType] = Success(WdlFileType)
  def read_lines(params: Seq[Try[WdlType]]): Try[WdlType] = Success(WdlArrayType(WdlStringType))
  def read_tsv(params: Seq[Try[WdlType]]): Try[WdlType] = Success(WdlArrayType(WdlArrayType(WdlStringType)))
  def read_map(params: Seq[Try[WdlType]]): Try[WdlType] = Success(WdlMapType(WdlStringType, WdlStringType))
  def read_object(params: Seq[Try[WdlType]]): Try[WdlType] = Success(WdlObjectType)
  def read_objects(params: Seq[Try[WdlType]]): Try[WdlType] = Success(WdlArrayType(WdlObjectType))
  def read_json(params: Seq[Try[WdlType]]): Try[WdlType] = Failure(new WdlExpressionException("Return type of read_json() can't be known statically"))
  def read_int(params: Seq[Try[WdlType]]): Try[WdlType] = Success(WdlIntegerType)
  def read_string(params: Seq[Try[WdlType]]): Try[WdlType] = Success(WdlStringType)
  def read_float(params: Seq[Try[WdlType]]): Try[WdlType] = Success(WdlFloatType)
  def read_boolean(params: Seq[Try[WdlType]]): Try[WdlType] = Success(WdlBooleanType)
  def write_lines(params: Seq[Try[WdlType]]): Try[WdlType] = Success(WdlFileType)
  def write_tsv(params: Seq[Try[WdlType]]): Try[WdlType] = Success(WdlFileType)
  def write_map(params: Seq[Try[WdlType]]): Try[WdlType] = Success(WdlFileType)
  def write_object(params: Seq[Try[WdlType]]): Try[WdlType] = Success(WdlFileType)
  def write_objects(params: Seq[Try[WdlType]]): Try[WdlType] = Success(WdlFileType)
  def write_json(params: Seq[Try[WdlType]]): Try[WdlType] = Success(WdlFileType)
  def glob(params: Seq[Try[WdlType]]): Try[WdlType] = Success(WdlArrayType(WdlFileType))
  def size(params: Seq[Try[WdlType]]): Try[WdlType] = Success(WdlFloatType)
  def sub(params: Seq[Try[WdlType]]): Try[WdlType] = Success(WdlStringType)
}

case object NoFunctions extends WdlStandardLibraryFunctions {
  override def glob(path: String, pattern: String): Seq[String] = throw new NotImplementedError()
  override def readFile(path: String): String = throw new NotImplementedError()
  override def writeTempFile(path: String, prefix: String, suffix: String, content: String): String = throw new NotImplementedError()
  override def stdout(params: Seq[Try[WdlValue]]): Try[WdlFile] = Failure(new NotImplementedError())
  override def stderr(params: Seq[Try[WdlValue]]): Try[WdlFile] = Failure(new NotImplementedError())
  override def read_json(params: Seq[Try[WdlValue]]): Try[WdlValue] = Failure(new NotImplementedError())
  override def write_tsv(params: Seq[Try[WdlValue]]): Try[WdlFile] = Failure(new NotImplementedError())
  override def write_json(params: Seq[Try[WdlValue]]): Try[WdlFile] = Failure(new NotImplementedError())
  override def size(params: Seq[Try[WdlValue]]): Try[WdlFloat] = Failure(new NotImplementedError())
  override def sub(params: Seq[Try[WdlValue]]): Try[WdlString] = Failure(new NotImplementedError())
}