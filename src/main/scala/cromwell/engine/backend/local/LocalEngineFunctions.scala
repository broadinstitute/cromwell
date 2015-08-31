package cromwell.engine.backend.local

import java.io.File
import java.nio.file.{Path, Paths}
import cromwell.binding.expression.WdlStandardLibraryFunctions
import cromwell.binding.types.{WdlArrayType, WdlFileType, WdlMapType, WdlStringType}
import cromwell.binding.values._
import cromwell.util.FileUtil
import cromwell.util.FileUtil.{EnhancedFile, EnhancedPath}

import scala.util.{Failure, Success, Try}

class LocalEngineFunctions(cwd: Path, stdout: Path, stderr: Path) extends WdlStandardLibraryFunctions {

  /**
   * Read the entire contents of a file from the specified `WdlValue`, where the file can be
   * specified either as a path via a `WdlString` (with magical handling of "stdout"), or
   * directly as a `WdlFile`.
   *
   * @throws UnsupportedOperationException for an unrecognized file reference, as this is intended
   *                                       to be wrapped in a `Try`.
   */
  private def fileContentsToString(value: WdlValue): String = {
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

  override protected def read_lines(params: Seq[Try[WdlValue]]): Try[WdlArray] = {
    for {
      singleArgument <- extractSingleArgument(params)
      lines = fileContentsToString(singleArgument).split("\n").map{WdlString}
    } yield WdlArray(WdlArrayType(WdlStringType), lines)
  }

  override protected def read_map(params: Seq[Try[WdlValue]]): Try[WdlMap] = {
    for {
      singleArgument <- extractSingleArgument(params)
      if singleArgument.wdlType == WdlFileType
      contents <- Success(Paths.get(singleArgument.asInstanceOf[WdlFile].valueString).slurp)
      wdlMap <- WdlMap.fromTsv(contents)
    } yield wdlMap
  }

  /**
   * Try to read an integer from the file referenced by the specified `WdlValue`.
   */
  override protected def read_int(params: Seq[Try[WdlValue]]): Try[WdlInteger] =
    read_string(params).map { s => WdlInteger(s.value.trim.toInt) }

  /**
   * Try to read a string from the file referenced by the specified `WdlValue`.
   */
  override protected def read_string(params: Seq[Try[WdlValue]]): Try[WdlString] = {
    for {
      singleArgument <- extractSingleArgument(params)
      string = fileContentsToString(singleArgument)
    } yield WdlString(string.stripSuffix("\n"))
  }

  override protected def write_lines(params: Seq[Try[WdlValue]]): Try[WdlFile] = {
    for {
      singleArgument <- extractSingleArgument(params)
      if singleArgument.wdlType.isInstanceOf[WdlArrayType]
      tsvSerialized <- singleArgument.asInstanceOf[WdlArray].tsvSerialize
      (path, writer) = FileUtil.tempFileAndWriter("array", cwd.toFile)
      _ <- Try(writer.write(tsvSerialized))
      _ <- Success(writer.close())
    } yield WdlFile(path.toAbsolutePath.toString)
  }

  override protected def write_map(params: Seq[Try[WdlValue]]): Try[WdlFile] = {
    for {
      singleArgument <- extractSingleArgument(params)
      if singleArgument.wdlType.isInstanceOf[WdlMapType]
      tsvSerialized <- singleArgument.asInstanceOf[WdlMap].tsvSerialize
      (path, writer) = FileUtil.tempFileAndWriter("map", cwd.toFile)
      _ <- Try(writer.write(tsvSerialized))
      _ <- Success(writer.close())
    } yield WdlFile(path.toAbsolutePath.toString)
  }
}
