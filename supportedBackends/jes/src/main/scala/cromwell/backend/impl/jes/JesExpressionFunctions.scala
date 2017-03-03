package cromwell.backend.impl.jes

import cromwell.backend.standard.{StandardExpressionFunctions, StandardExpressionFunctionsParams}
import cromwell.filesystems.gcs.GcsPathBuilder

class JesExpressionFunctions(standardParams: StandardExpressionFunctionsParams)
  extends StandardExpressionFunctions(standardParams) {

  override def preMapping(str: String) =
    if (!GcsPathBuilder.isValidGcsUrl(str)) {
      callContext.root.resolve(str.stripPrefix("/")).pathAsString
    } else {
      str
    }
}
