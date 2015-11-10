package cromwell.engine.backend.jes

import cromwell.binding.expression.WdlStandardLibraryFunctions
import cromwell.binding.types.{WdlArrayType, WdlFileType, WdlObjectType, WdlStringType}
import cromwell.binding.values._
import cromwell.engine.backend.jes.authentication.ProductionJesAuthentication
import cromwell.util.google.GoogleCloudStoragePath

import scala.language.postfixOps
import scala.util.{Success, Try}

/**
 * Implementation of WDL standard library functions for the JES backend.
 */
class JesEngineFunctionsWithoutCallContext extends WdlStandardLibraryFunctions with ProductionJesAuthentication {

  private def readFromPath(value: String): String = {
    // .get here because engine functions should throw exception if they fail.  Evaluator will catch it
    authenticated { _.storage.slurpFile(GoogleCloudStoragePath.parse(value).get) }
  }

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
      case f: WdlFile => readFromPath(f.value)
      case f: WdlString => readFromPath(f.value)
      case e => throw new UnsupportedOperationException("Unsupported argument " + e + " (expected JES URI)")
    }
  }
}

/**
 * Implementation of WDL standard library functions for the JES backend.
 */
class JesEngineFunctions(jesBackendCall: JesBackendCall) extends JesEngineFunctionsWithoutCallContext {

  override protected def glob(params: Seq[Try[WdlValue]]): Try[WdlArray] = {
    for {
      singleArgument <- extractSingleArgument(params)
      files = authenticated { _.storage.listContents(jesBackendCall.globOutputPath(singleArgument.valueString)) }
      wdlFiles = files map { WdlFile(_, isGlob = false) }
    } yield WdlArray(WdlArrayType(WdlFileType), wdlFiles toSeq)
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
