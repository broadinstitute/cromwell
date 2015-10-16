package cromwell.engine.backend.local

import java.io.File
import java.nio.file.{Path, Paths}

import cromwell.binding.expression.WdlStandardLibraryFunctions
import cromwell.binding.types.WdlArrayType._
import cromwell.binding.types._
import cromwell.binding.values._
import cromwell.util.FileUtil
import cromwell.util.FileUtil.{EnhancedFile, EnhancedPath}

import scala.util.{Failure, Success, Try}

class LocalEngineFunctionsWithoutCallContext extends WdlStandardLibraryFunctions {
  protected def fileContentsToString(value: WdlValue): String = {
    value match {
      case f: WdlFile => new File(f.value).slurp
      case s: WdlString => Paths.get(s.value).slurp
      case e => throw new UnsupportedOperationException("Unsupported argument " + e)
    }
  }

  override protected def read_lines(params: Seq[Try[WdlValue]]): Try[WdlArray] = {
    for {
      singleArgument <- extractSingleArgument(params)
      lines = fileContentsToString(singleArgument).split("\n").map{WdlString}
    } yield WdlArray(WdlArrayType(WdlStringType), lines)
  }

  override protected def read_map(params: Seq[Try[WdlValue]]): Try[WdlMap] = {
    for {
      singleArgument <- extractSingleArgument(params)
      contents <- Success(fileContentsToString(singleArgument))
      wdlMap <- WdlMap.fromTsv(contents)
    } yield wdlMap
  }

  private def extractObjectArray(params: Seq[Try[WdlValue]]): Try[Array[WdlObject]] = for {
    singleArgument <- extractSingleArgument(params)
    contents <- Success(fileContentsToString(singleArgument))
    wdlObjects <- WdlObject.fromTsv(contents)
  } yield wdlObjects

  override protected def read_object(params: Seq[Try[WdlValue]]): Try[WdlObject] = {
    extractObjectArray(params) map {
      case array if array.length == 1 => array.head
      case _ => throw new IllegalArgumentException("read_object yields an Object and thus can only read 2-rows TSV files. Try using read_objects instead.")
    }
  }

  override def read_objects(params: Seq[Try[WdlValue]]): Try[WdlArray] = {
    extractObjectArray(params) map { WdlArray(WdlArrayType(WdlObjectType), _) }
  }

  /**
   * Try to read a string from the file referenced by the specified `WdlValue`.
   */
  override protected def read_string(params: Seq[Try[WdlValue]]): Try[WdlString] = {
    for {
      singleArgument <- extractSingleArgument(params)
      string = fileContentsToString(singleArgument)
    } yield WdlString(string.stripSuffix("\n"))
  }
}

class LocalEngineFunctions(cwd: Path, stdout: Path, stderr: Path) extends LocalEngineFunctionsWithoutCallContext {

  /**
   * Read the entire contents of a file from the specified `WdlValue`, where the file can be
   * specified either as a path via a `WdlString` (with magical handling of "stdout"), or
   * directly as a `WdlFile`.
   *
   * @throws UnsupportedOperationException for an unrecognized file reference, as this is intended
   *                                       to be wrapped in a `Try`.
   */
  override def fileContentsToString(value: WdlValue): String = {
    value match {
      case f: WdlFile => new File(f.value).slurp
      case s: WdlString => cwd.resolve(s.value).slurp
      case e => throw new UnsupportedOperationException("Unsupported argument " + e)
    }
  }

  override protected def stdout(params: Seq[Try[WdlValue]]): Try[WdlFile] = {
    if (params.nonEmpty) {
      Failure(new UnsupportedOperationException("stdout() takes zero parameters"))
    } else {
      Success(WdlFile(stdout.toAbsolutePath.toString))
    }
  }

  override protected def stderr(params: Seq[Try[WdlValue]]): Try[WdlFile] = {
    if (params.nonEmpty) {
      Failure(new UnsupportedOperationException("stderr() takes zero parameters"))
    } else {
      Success(WdlFile(stderr.toAbsolutePath.toString))
    }
  }

  override protected def write_lines(params: Seq[Try[WdlValue]]): Try[WdlFile] = {
    for {
      singleArgument <- extractSingleArgument(params)
      if singleArgument.wdlType.isInstanceOf[WdlArrayType]
      tsvSerialized <- singleArgument.asInstanceOf[WdlArray].tsvSerialize
      file <- writeContent("array", tsvSerialized)
    } yield file
  }

  override protected def write_map(params: Seq[Try[WdlValue]]): Try[WdlFile] = {
    for {
      singleArgument <- extractSingleArgument(params)
      if singleArgument.wdlType.isInstanceOf[WdlMapType]
      tsvSerialized <- singleArgument.asInstanceOf[WdlMap].tsvSerialize
      file <- writeContent("map", tsvSerialized)
    } yield file
  }

  override protected def write_object(params: Seq[Try[WdlValue]]): Try[WdlFile] = {
    for {
      singleArgument <- extractSingleArgument(params)
      if singleArgument.wdlType == WdlObjectType
      tsvSerialized <- singleArgument.asInstanceOf[WdlObject].tsvSerialize
      file <- writeContent("object", tsvSerialized)
    } yield file
  }

  override protected def write_objects(params: Seq[Try[WdlValue]]): Try[WdlFile] = {
    for {
      singleArgument <- extractSingleArgument(params)
      if singleArgument.wdlType.isAnArrayOf(WdlObjectType)
      tsvSerialized <- singleArgument.asInstanceOf[WdlArray].tsvSerialize
      file <- writeContent("array", tsvSerialized)
    } yield file
  }

  protected def writeContent(baseName: String, content: String) = {
    val (path, writer) = FileUtil.tempFileAndWriter("array", cwd.toFile)
    try {
      writer.write(content)
      Success(WdlFile(path.toAbsolutePath.toString))
    } catch {
      case t: Throwable => Failure(t)
    } finally {
      writer.close()
    }
  }
}
