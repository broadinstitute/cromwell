package cromwell.engine.backend.jes

import cromwell.binding.expression.WdlStandardLibraryFunctions
import cromwell.binding.types.{WdlArrayType, WdlFileType, WdlObjectType, WdlStringType}
import cromwell.binding.values._
import cromwell.engine.backend.jes.authentication.ProductionJesAuthentication
import cromwell.util.google.GoogleCloudStoragePath

import scala.language.postfixOps
import scala.util.{Success, Try}

object JesEngineFunctions {
  implicit class UriString(val str: String) extends AnyVal {
    def isGcsUrl: Boolean = str.startsWith("gs://")
    def isUriWithProtocol: Boolean = "^[a-zA-Z]+://".r.findFirstIn(str).nonEmpty
  }
}

/**
 * Implementation of WDL standard library functions for the JES backend.
 */
class JesEngineFunctionsWithoutCallContext extends WdlStandardLibraryFunctions with ProductionJesAuthentication {

  import JesEngineFunctions._

  protected def readFromGcsUri(value: String): String = {
    // .get here because engine functions should throw exception if they fail.  Evaluator will catch it
    authenticated { _.storage.slurpFile(GoogleCloudStoragePath.parse(value).get) }
  }

  /**
    * Read the entire contents of a file from the specified `WdlValue`.  the `WdlValue` must be
    * either a `WdlString` or `WdlFile` (i.e. `WdlStringLike`) and must be a full gs:// url
    *
    * @throws UnsupportedOperationException for an unrecognized file reference, as this is intended
    *                                       to be wrapped in a `Try`.
    */
  override def fileContentsToString(path: String): String = path match {
    case s if s.isGcsUrl => readFromGcsUri(s)
    case s if s.isUriWithProtocol => throw new UnsupportedOperationException(s"URI Scheme not supported: $s")
    case e => throw new UnsupportedOperationException("Unsupported argument " + e + " (expected GCS URI)")
  }
}

/**
 * Implementation of WDL standard library functions for the JES backend.
 */
class JesEngineFunctions(jesBackendCall: JesBackendCall) extends JesEngineFunctionsWithoutCallContext {

  import JesEngineFunctions._

  /**
    * Read the entire contents of a file from the specified `WdlValue`.  the `WdlValue` must be
    * either a `WdlString` or `WdlFile` (i.e. `WdlStringLike`) and must be a full gs:// url OR
    * a relative path in the call output directory in which the call's GCS path will be prepended
    *
    * @throws UnsupportedOperationException for an unrecognized file reference, as this is intended
    *                                       to be wrapped in a `Try`.
    */
  override def fileContentsToString(path: String): String = path match {
    case s if s.isGcsUrl => readFromGcsUri(s)
    case s if s.isUriWithProtocol => throw new UnsupportedOperationException(s"URI Scheme not supported: $s")
    case s => readFromGcsUri(s"${jesBackendCall.callGcsPath}/$s")
  }


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
