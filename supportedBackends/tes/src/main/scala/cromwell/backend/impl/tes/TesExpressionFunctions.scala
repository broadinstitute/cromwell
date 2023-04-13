package cromwell.backend.impl.tes

import cromwell.backend.sfs.SharedFileSystemExpressionFunctions
import cromwell.backend.standard.StandardExpressionFunctionsParams

class TesExpressionFunctions(standardParams: StandardExpressionFunctionsParams) extends SharedFileSystemExpressionFunctions(standardParams) {
  override lazy implicit val ec = standardParams.executionContext;

  override def preMapping(str: String) = {
    if (str.startsWith("/") || str.startsWith("ftp://")) str
    else standardParams.callContext.root.resolve(str).toString
  }
}
