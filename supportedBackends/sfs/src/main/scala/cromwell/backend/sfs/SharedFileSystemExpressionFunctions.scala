package cromwell.backend.sfs

import cromwell.backend.io._
import cromwell.backend.standard.{DefaultStandardExpressionFunctionsParams, StandardExpressionFunctions, StandardExpressionFunctionsParams}
import cromwell.core.CallContext
import cromwell.core.path.{DefaultPath, Path, PathBuilder}

object SharedFileSystemExpressionFunctions {
  def apply(jobPaths: JobPaths, pathBuilders: List[PathBuilder]): SharedFileSystemExpressionFunctions = {
    new SharedFileSystemExpressionFunctions(pathBuilders, jobPaths.callContext)
  }
}

class SharedFileSystemExpressionFunctions(standardParams: StandardExpressionFunctionsParams)
  extends StandardExpressionFunctions(standardParams) {

  def this(pathBuilders: List[PathBuilder], callContext: CallContext) = {
    this(DefaultStandardExpressionFunctionsParams(pathBuilders, callContext))
  }

  override def postMapping(path: Path) = {
    path match {
      case _: DefaultPath if !path.isAbsolute => callContext.root.resolve(path)
      case _ => path
    }
  }
}
