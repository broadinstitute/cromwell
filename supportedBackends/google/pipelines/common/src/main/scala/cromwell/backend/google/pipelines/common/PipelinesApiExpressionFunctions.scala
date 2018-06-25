package cromwell.backend.google.pipelines.common

import cromwell.backend.standard.{StandardExpressionFunctions, StandardExpressionFunctionsParams}
import cromwell.core.CallContext
import cromwell.core.io.{CallCorePathFunctionSet, IoCommandBuilder}
import cromwell.core.path.PathFactory.PathBuilders
import cromwell.filesystems.gcs.GcsPathBuilder
import cromwell.filesystems.gcs.GcsPathBuilder.{InvalidGcsPath, PossiblyValidRelativeGcsPath, ValidFullGcsPath}
import cromwell.filesystems.gcs.batch.GcsBatchCommandBuilder

class PipelinesApiPathFunctions(pathBuilders: PathBuilders, callContext: CallContext) extends CallCorePathFunctionSet(pathBuilders, callContext) {
  override def relativeToHostCallRoot(path: String) = {
    GcsPathBuilder.validateGcsPath(path) match {
      case _: ValidFullGcsPath => path
      case _ => callContext.root.resolve(path.stripPrefix("file://").stripPrefix("/")).pathAsString
    }
  }
}

class PipelinesApiExpressionFunctions(standardParams: StandardExpressionFunctionsParams)
  extends StandardExpressionFunctions(standardParams) {
  override lazy val ioCommandBuilder: IoCommandBuilder = GcsBatchCommandBuilder

  override def preMapping(str: String) = {
    GcsPathBuilder.validateGcsPath(str) match {
      case _: ValidFullGcsPath => str
      case PossiblyValidRelativeGcsPath => callContext.root.resolve(str.stripPrefix("/")).pathAsString
      case invalid: InvalidGcsPath => throw new IllegalArgumentException(invalid.errorMessage)
    }
  }

  override lazy val pathFunctions = new PipelinesApiPathFunctions(pathBuilders, callContext)
}
