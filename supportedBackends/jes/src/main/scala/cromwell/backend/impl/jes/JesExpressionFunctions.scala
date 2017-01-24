package cromwell.backend.impl.jes

import cromwell.backend.standard.{StandardExpressionFunctions, StandardExpressionFunctionsParams}
import cromwell.core.path.PathImplicits._
import cromwell.filesystems.gcs.GcsPathBuilder

class JesExpressionFunctions(standardParams: StandardExpressionFunctionsParams)
  extends StandardExpressionFunctions(standardParams) {

  override def preMapping(str: String): String =
    if (!GcsPathBuilder.isValidGcsUrl(str)) {
      callContext.root.resolve(str.stripPrefix("/")).toRealString
    } else {
      str
    }
}
