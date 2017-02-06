package cromwell.backend.sfs

import java.nio.file.Path

import cromwell.backend.io._
import cromwell.backend.standard.{DefaultStandardExpressionFunctionsParams, StandardExpressionFunctions, StandardExpressionFunctionsParams}
import cromwell.core.CallContext
import cromwell.core.path.PathBuilder

object SharedFileSystemExpressionFunctions {
  private val LocalFileSystemScheme = "file"

  def isLocalPath(path: Path): Boolean =
    path.toUri.getScheme == SharedFileSystemExpressionFunctions.LocalFileSystemScheme

  def apply(jobPaths: JobPaths, pathBuilders: List[PathBuilder]): SharedFileSystemExpressionFunctions = {
    new SharedFileSystemExpressionFunctions(pathBuilders, jobPaths.callContext)
  }
}

class SharedFileSystemExpressionFunctions(standardParams: StandardExpressionFunctionsParams)
  extends StandardExpressionFunctions(standardParams) {

  def this(pathBuilders: List[PathBuilder], callContext: CallContext) = {
    this(DefaultStandardExpressionFunctionsParams(pathBuilders, callContext))
  }

  override def postMapping(path: Path): Path =
    if (!path.isAbsolute && SharedFileSystemExpressionFunctions.isLocalPath(path))
      callContext.root.resolve(path)
    else
      path
}
