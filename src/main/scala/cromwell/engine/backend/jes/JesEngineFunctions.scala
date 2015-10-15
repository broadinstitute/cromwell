package cromwell.engine.backend.jes

import cromwell.binding.expression.WdlStandardLibraryFunctions
import cromwell.binding.types.{WdlObjectType, WdlArrayType, WdlStringType}
import cromwell.binding.values._
import cromwell.util.google.{GoogleCloudStorage, GoogleCloudStoragePath}

import scala.util.{Failure, Success, Try}

class JesEngineFunctionsWithoutCallContext(storage: GoogleCloudStorage) extends WdlStandardLibraryFunctions {

  protected def readFromPath(value: String): String = {
    // .get here because engine functions should throw exception if they fail.  Evaluator will catch it
    storage.slurpFile(GoogleCloudStoragePath.parse(value).get)
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

  /**
   * Read all lines from the file referenced by the first parameter
   */
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

  private def extractObjectOrArray(params: Seq[Try[WdlValue]]) = for {
    singleArgument <- extractSingleArgument(params)
    contents <- Success(fileContentsToString(singleArgument))
    wdlObjects <- WdlObject.fromTsv(contents)
  } yield wdlObjects

  override protected def read_object(params: Seq[Try[WdlValue]]): Try[WdlObject] = {
    extractObjectOrArray(params) map {
      case array if array.length == 1 => array.head
      case _ => throw new IllegalArgumentException("read_object yields an Object and thus can only read 2-rows TSV files. Try using read_objects instead.")
    }
  }

  override def read_objects(params: Seq[Try[WdlValue]]): Try[WdlArray] = {
    extractObjectOrArray(params) map { WdlArray(WdlArrayType(WdlObjectType), _) }
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

/**
 * Implementation of WDL standard library functions for the JES backend.
 */
class JesEngineFunctions(jesBackendCall: JesBackendCall) extends JesEngineFunctionsWithoutCallContext(jesBackendCall.jesConnection.storage) {

  override def readFromPath(value: String): String = {
    val gcsPath = GoogleCloudStoragePath.parse(value) match {
      case Success(path) => path
      case Failure(ex) => GoogleCloudStoragePath(jesBackendCall.callDir + s"/$value")
    }
    jesBackendCall.jesConnection.storage.slurpFile(gcsPath)
  }

  override protected def stdout(params: Seq[Try[WdlValue]]): Try[WdlFile] = {
    val newPath = GoogleCloudStoragePath(jesBackendCall.stdoutJesOutput.gcs)
    Success(WdlFile(newPath.toString))
  }

  override protected def stderr(params: Seq[Try[WdlValue]]): Try[WdlFile] = {
    val newPath = GoogleCloudStoragePath(jesBackendCall.stderrJesOutput.gcs)
    Success(WdlFile(newPath.toString))
  }
}
