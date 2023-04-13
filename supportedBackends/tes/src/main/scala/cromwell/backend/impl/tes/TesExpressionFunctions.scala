package cromwell.backend.impl.tes

import cromwell.backend.sfs.SharedFileSystemExpressionFunctions
import cromwell.backend.standard.StandardExpressionFunctionsParams

import scala.concurrent.Future

class TesExpressionFunctions(standardParams: StandardExpressionFunctionsParams) extends SharedFileSystemExpressionFunctions(standardParams) {
  override lazy implicit val ec = standardParams.executionContext;

  override def preMapping(str: String) = {
    if (str.startsWith("/") || str.startsWith("ftp://")) str
    else standardParams.callContext.root.resolve(str).toString
  }

  override def glob(pattern: String): Future[Seq[String]] = {
    import wom.values.GlobFunctions._
    val globPatternName = globName(pattern)
    val listFilePath = callContext.root.resolve(s"${globName(pattern)}.list")
    asyncIo.readLinesAsync(listFilePath.toAbsolutePath) map { lines =>
      lines.toList map { fileName =>
        (callContext.root /  globPatternName  / fileName).pathAsString
      }
    }
  }
}
