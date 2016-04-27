package cromwell.backend

import java.nio.file.Path

import better.files._
import cromwell.backend.validation.MemorySize
import wdl4s.expression.WdlStandardLibraryFunctions
import wdl4s.parser.MemoryUnit
import wdl4s.types._
import wdl4s.values._
import wdl4s.TsvSerializable

import scala.language.postfixOps
import scala.util.{Failure, Success, Try}
import scala.language.implicitConversions

trait WdlStandardLibraryImpl extends WdlStandardLibraryFunctions {

  def toPath(str: String): Path

  override final def fileContentsToString(path: String): String = toPath(path).toAbsolutePath.contentAsString

  private def writeContent(baseName: String, content: String): Try[WdlFile] = {
    Try(WdlFile(writeTempFile(tempFilePath, s"$baseName.", ".tmp", content)))
  }

  /**
    * Asserts that the parameter list contains a single parameter which will be interpreted
    * as a File and attempts to read the contents of that file
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

  override def readFile(path: String): String = toPath(path).toAbsolutePath.contentAsString

  /**
    * Read all lines from the file referenced by the first parameter and return an Array[String]
    */
  override def read_lines(params: Seq[Try[WdlValue]]): Try[WdlArray] = {
    for {
      contents <- readContentsFromSingleFileParameter(params)
      lines = contents.split("\n")
    } yield WdlArray(WdlArrayType(WdlStringType), lines map WdlString)
  }

  override def read_map(params: Seq[Try[WdlValue]]): Try[WdlMap] = {
    for {
      contents <- readContentsFromSingleFileParameter(params)
      wdlMap <- WdlMap.fromTsv(contents)
    } yield wdlMap
  }

  override def read_object(params: Seq[Try[WdlValue]]): Try[WdlObject] = {
    extractObjects(params) map {
      case array if array.length == 1 => array.head
      case _ => throw new IllegalArgumentException("read_object yields an Object and thus can only read 2-rows TSV files. Try using read_objects instead.")
    }
  }

  override def read_objects(params: Seq[Try[WdlValue]]): Try[WdlArray] = extractObjects(params) map { WdlArray(WdlArrayType(WdlObjectType), _) }

  /**
    * Try to read a string from the file referenced by the specified `WdlValue`.
    */
  override def read_string(params: Seq[Try[WdlValue]]): Try[WdlString] = readContentsFromSingleFileParameter(params).map(s => WdlString(s.trim))

  /**
    * Read a file in TSV format into an Array[Array[String]]
    */
  override def read_tsv(params: Seq[Try[WdlValue]]): Try[WdlArray] = {
    for {
      contents <- readContentsFromSingleFileParameter(params)
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

  override def sub(params: Seq[Try[WdlValue]]): Try[WdlString] = {
    def extractArguments = params.size match {
      case 3 => Success((params.head, params(1), params(2)))
      case n => Failure(new IllegalArgumentException(s"Invalid number of parameters for engine function sub: $n. sub takes exactly 3 parameters."))
    }

    def validateArguments(values: (Try[WdlValue], Try[WdlValue], Try[WdlValue])) = values match {
      case (Success(strValue), Success(WdlString(pattern)), Success(replaceValue))
        if WdlStringType.isCoerceableFrom(strValue.wdlType) &&
           WdlStringType.isCoerceableFrom(replaceValue.wdlType) =>
        Success((strValue.valueString, pattern, replaceValue.valueString))
      case _ => Failure(new IllegalArgumentException(s"Invalid parameters for engine function sub: $values."))
    }

    for {
      args <- extractArguments
      (str, pattern, replace) <- validateArguments(args)
    } yield WdlString(pattern.r.replaceAllIn(str, replace))
  }

  /**
    * Try to read a boolean from the file referenced by the specified `WdlValue`.
    */
  override def read_boolean(params: Seq[Try[WdlValue]]): Try[WdlBoolean] =
    read_string(params) map { s => WdlBoolean(java.lang.Boolean.parseBoolean(s.value.trim.toLowerCase)) }

  private def writeToTsv(params: Seq[Try[WdlValue]], wdlClass: Class[_ <: WdlValue with TsvSerializable]) = {
    for {
      singleArgument <- extractSingleArgument(params)
      downcast <- Try(wdlClass.cast(singleArgument))
      tsvSerialized <- downcast.tsvSerialize
      file <- writeContent(wdlClass.getSimpleName.toLowerCase, tsvSerialized)
    } yield file
  }

  override def size(params: Seq[Try[WdlValue]]): Try[WdlFloat] = {
    def toUnit(wdlValue: Try[WdlValue]) = wdlValue flatMap { unit => Try(MemoryUnit.fromSuffix(unit.valueString)) }

    def fileSize(wdlValue: Try[WdlValue], convertTo: Try[MemoryUnit] = Success(MemoryUnit.Bytes)) = {
      for {
        value <- wdlValue
        unit <- convertTo
      } yield MemorySize(toPath(value.valueString).toAbsolutePath.size.toDouble, MemoryUnit.Bytes).to(unit).amount
    }

    params match {
      case oneArg if params.length == 1 => fileSize(params.head) map WdlFloat.apply
      case twoArgs if params.length == 2 => fileSize(params.head, toUnit(params(1))) map WdlFloat.apply
      case _ => Failure(new UnsupportedOperationException(s"Expected one or two parameters but got ${params.length} instead."))
    }
  }

  override def write_lines(params: Seq[Try[WdlValue]]): Try[WdlFile] = writeToTsv(params, classOf[WdlArray])
  override def write_map(params: Seq[Try[WdlValue]]): Try[WdlFile] = writeToTsv(params, classOf[WdlMap])
  override def write_object(params: Seq[Try[WdlValue]]): Try[WdlFile] = writeToTsv(params, classOf[WdlObject])
  override def write_objects(params: Seq[Try[WdlValue]]): Try[WdlFile] = writeToTsv(params, classOf[WdlArray])
  override def write_tsv(params: Seq[Try[WdlValue]]): Try[WdlFile] = writeToTsv(params, classOf[WdlArray])
  override def writeTempFile(path: String, prefix: String, suffix: String, content: String): String = {
    // This may be called multiple times with the same inputs.  Calling this twice with the same
    // parameters should yield the same return value.
    val rootDir = toPath(path).toAbsolutePath
    rootDir.createDirectories

    val fullPath = rootDir.resolve(s"$prefix${content.md5Sum}$suffix")
    if (!fullPath.exists) fullPath.write(content)

    fullPath.toString
  }

  override def glob(params: Seq[Try[WdlValue]]): Try[WdlArray] = {
    for {
      singleArgument <- extractSingleArgument(params)
      globVal = singleArgument.valueString
      files = glob(globPath(globVal), globVal)
      wdlFiles = files map { WdlFile(_, isGlob = false) }
    } yield WdlArray(WdlArrayType(WdlFileType), wdlFiles toSeq)
  }

  private def fail(name: String) = Failure(new NotImplementedError(s"$name() not implemented yet"))
  override def read_json(params: Seq[Try[WdlValue]]): Try[WdlValue] = fail("read_json")
  override def write_json(params: Seq[Try[WdlValue]]): Try[WdlFile] = fail("write_json")
  override def stdout(params: Seq[Try[WdlValue]]): Try[WdlFile] = fail("stdout")
  override def stderr(params: Seq[Try[WdlValue]]): Try[WdlFile] = fail("stderr")
}
