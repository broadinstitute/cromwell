package cromwell.engine.backend.local

import java.nio.file.Path

import better.files._
import cromwell.binding.IOInterface
import cromwell.binding.expression.WdlStandardLibraryFunctions
import cromwell.binding.types.WdlArrayType._
import cromwell.binding.types._
import cromwell.binding.values._

import scala.language.postfixOps
import scala.util.{Failure, Success, Try}

class LocalEngineFunctionsWithoutCallContext(interface: IOInterface) extends WdlStandardLibraryFunctions {
  override def fileContentsToString(path: String): String = interface.readFile(path)
}

class LocalEngineFunctions(cwd: Path, stdout: Path, stderr: Path, interface: IOInterface) extends LocalEngineFunctionsWithoutCallContext(interface) {

  /**
   * Read the entire contents of a file from the specified `WdlValue`, where the file can be
   * specified either as a path via a `WdlString` (with magical handling of "stdout"), or
   * directly as a `WdlFile`.
   *
   * @throws UnsupportedOperationException for an unrecognized file reference, as this is intended
   *                                       to be wrapped in a `Try`.
   */
  override def fileContentsToString(path: String): String = interface.readFile(cwd.resolve(path).toString)

  override protected def stdout(params: Seq[Try[WdlValue]]): Try[WdlFile] = {
    if (params.nonEmpty) {
      Failure(new UnsupportedOperationException("stdout() takes zero parameters"))
    } else {
      Success(WdlFile(stdout.fullPath))
    }
  }

  override protected def stderr(params: Seq[Try[WdlValue]]): Try[WdlFile] = {
    if (params.nonEmpty) {
      Failure(new UnsupportedOperationException("stderr() takes zero parameters"))
    } else {
      Success(WdlFile(stderr.fullPath))
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

  override protected def glob(params: Seq[Try[WdlValue]]): Try[WdlArray] = {
    for {
      singleArgument <- extractSingleArgument(params)
      if singleArgument.wdlType == WdlFileType || singleArgument.wdlType == WdlStringType
      files <- filesMatchingGlob(singleArgument.valueString)
    } yield WdlArray(WdlArrayType(WdlFileType), files)
  }

  private def filesMatchingGlob(glob: String): Try[Seq[WdlValue]] = Try {
    interface.glob(("." / cwd.toString).fullPath, s"**/$glob") map { WdlFile(_) }
  }

  protected def writeContent(baseName: String, content: String): Try[WdlFile] = {
    Try(WdlFile(interface.writeTempFile(cwd.toString, s"$baseName.", ".tmp", content)))
  }
}
