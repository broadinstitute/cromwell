package cromwell.engine.backend.jes

import cromwell.binding.values._
import cromwell.engine.EngineFunctions
import cromwell.util.GcsUtil

import scala.util.{Success, Try}

/**
 * Implementation of engine functions for the JES backend.
 */
class JesEngineFunctions(secretsFile: String) extends EngineFunctions {

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
      case f: WdlGcsObject => GcsUtil.slurp(f.value, secretsFile)
      case e => throw new UnsupportedOperationException("Unsupported argument " + e + " (expected JES URI)")
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

  //TODO: Implement correctly
  override protected def stdout(params: Seq[Try[WdlValue]]): Try[WdlGcsObject] = {
    Success(WdlGcsObject("gs://chrisl-dsde-dev/stdout"))
  }

  //TODO: Implement correctly
  override protected def stderr(params: Seq[Try[WdlValue]]): Try[WdlGcsObject] = {
    Success(WdlGcsObject("gs://chrisl-dsde-dev/stderr"))
  }
}
