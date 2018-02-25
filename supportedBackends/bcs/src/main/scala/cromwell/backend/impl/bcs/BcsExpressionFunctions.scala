package cromwell.backend.impl.bcs

import cromwell.backend.standard.{StandardExpressionFunctions, StandardExpressionFunctionsParams}
import cromwell.filesystems.oss.OssPathBuilder
import cromwell.filesystems.oss.OssPathBuilder.{InvalidOssPath, PossiblyValidRelativeOssPath, ValidFullOssPath}

final case class BcsExpressionFunctions(override val standardParams: StandardExpressionFunctionsParams)
  extends StandardExpressionFunctions(standardParams) {

  override def preMapping(str: String) = {
    OssPathBuilder.validateOssPath(str) match {
      case _: ValidFullOssPath => str
      case PossiblyValidRelativeOssPath => callContext.root.resolve(str.stripPrefix("/")).pathAsString
      case invalid: InvalidOssPath => throw new IllegalArgumentException(invalid.errorMessage)
    }
  }
}