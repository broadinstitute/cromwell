package cromwell.engine.backend.local

import java.nio.file.Paths

import cromwell.binding.values.{WdlFile, WdlInteger, WdlString, WdlValue}
import cromwell.engine.EngineFunctions
import cromwell.util.FileUtil

import scala.util.{Success, Failure, Try}


class LocalEngineFunctions(executionContext: TaskExecutionContext) extends EngineFunctions {

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
      case f: WdlFile => FileUtil.slurp(Paths.get(f.value))
      case e => throw new UnsupportedOperationException("Unsupported argument " + e)
    }
  }

  /**
   * Try to read a string from the file referenced by the specified `WdlValue`.
   */
  override protected def read_string(params: Seq[Try[WdlValue]]): Try[WdlString] = {
    for {
      singleArgument <- extractSingleArgument(params)
      string = fileContentsToString(singleArgument)
    } yield WdlString(string)
  }

  /**
   * Try to read an integer from the file referenced by the specified `WdlValue`.
   */
  override protected def read_int(params: Seq[Try[WdlValue]]): Try[WdlInteger] =
    read_string(params).map { s => WdlInteger(s.value.trim.toInt) }

  override protected def stdout(params: Seq[Try[WdlValue]]): Try[WdlFile] = {
    if (params.nonEmpty) {
      Failure(new UnsupportedOperationException("stdout() takes zero parameters"))
    } else {
      Success(WdlFile(executionContext.stdout.toAbsolutePath.toString))
    }
  }

  override protected def stderr(params: Seq[Try[WdlValue]]): Try[WdlFile] = {
    if (params.nonEmpty) {
      Failure(new UnsupportedOperationException("stderr() takes zero parameters"))
    } else {
      Success(WdlFile(executionContext.stderr.toAbsolutePath.toString))
    }
  }
}
