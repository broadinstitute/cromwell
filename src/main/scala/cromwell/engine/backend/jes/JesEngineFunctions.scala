package cromwell.engine.backend.jes

import cromwell.binding.expression.WdlStandardLibraryFunctions
import cromwell.binding.types.{WdlArrayType, WdlStringType}
import cromwell.binding.values._
import cromwell.util.google.GoogleCloudStoragePath

import scala.util.{Failure, Success, Try}

/**
 * Implementation of WDL standard library functions for the JES backend.
 */
case class JesEngineFunctions(jesBackendCall: JesBackendCall) extends WdlStandardLibraryFunctions {

  private def readFromPath(value: String): String = {
    val gcsPath = GoogleCloudStoragePath.parse(value) match {
      case Success(path) => path
      case Failure(ex) => GoogleCloudStoragePath(jesBackendCall.callDir + s"/$value")
    }
    jesBackendCall.jesConnection.storage.slurpFile(gcsPath)
  }
  
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
      case f: WdlFile => readFromPath(f.value)
      case f: WdlString => readFromPath(f.value)
      case e => throw new UnsupportedOperationException("Unsupported argument " + e + " (expected JES URI)")
    }
  }

  override protected def stdout(params: Seq[Try[WdlValue]]): Try[WdlFile] = {
    val newPath = GoogleCloudStoragePath(jesBackendCall.callDir.bucket, jesBackendCall.callDir.objectName + "/" + JesBackend.LocalStdoutValue)
    Success(WdlFile(newPath.toString))
  }

  override protected def stderr(params: Seq[Try[WdlValue]]): Try[WdlFile] = {
    val newPath = GoogleCloudStoragePath(jesBackendCall.callDir.bucket, jesBackendCall.callDir.objectName + "/" + JesBackend.LocalStderrValue)
    Success(WdlFile(newPath.toString))
  }

  override protected def read_lines(params: Seq[Try[WdlValue]]): Try[WdlArray] = {
    for {
      singleArgument <- extractSingleArgument(params)
      lines = fileContentsToString(singleArgument).split("\n").map{WdlString}
    } yield WdlArray(WdlArrayType(WdlStringType), lines)
  }

  /**
   * Try to read an integer from the file referenced by the specified `WdlValue`.
   */
  override protected def read_int(params: Seq[Try[WdlValue]]): Try[WdlInteger] = {
    read_string(params) map { s => WdlInteger(s.value.trim.toInt) }
  }

  /**
   * Try to read an float from the file referenced by the specified `WdlValue`.
   */
  override protected def read_float(params: Seq[Try[WdlValue]]): Try[WdlFloat] = {
    read_string(params) map { s => WdlFloat(s.value.trim.toDouble) }
  }

  /**
   * Try to read a string from the file referenced by the specified `WdlValue`.
   */
  override protected def read_string(params: Seq[Try[WdlValue]]): Try[WdlString] = {
    for {
      singleArgument <- extractSingleArgument(params)
      string = fileContentsToString(singleArgument)
    } yield WdlString(string.trim)
  }
}
