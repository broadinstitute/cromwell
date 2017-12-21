package cromwell.backend.sfs

import akka.actor.ActorRef
import cromwell.backend.io._
import cromwell.backend.standard.{DefaultStandardExpressionFunctionsParams, StandardExpressionFunctions, StandardExpressionFunctionsParams}
import cromwell.core.CallContext
import cromwell.core.path.{DefaultPath, Path, PathBuilder}

import scala.concurrent.ExecutionContext

object SharedFileSystemExpressionFunctions {
  def apply(jobPaths: JobPaths,
            pathBuilders: List[PathBuilder],
            ioActorEndpoint: ActorRef,
            ec: ExecutionContext): SharedFileSystemExpressionFunctions = {
    new SharedFileSystemExpressionFunctions(pathBuilders, jobPaths.callContext, ioActorEndpoint, ec)
  }
}

class SharedFileSystemExpressionFunctions(standardParams: StandardExpressionFunctionsParams)
  extends StandardExpressionFunctions(standardParams) {

  def this(pathBuilders: List[PathBuilder],
           callContext: CallContext,
           ioActorEndpoint: ActorRef,
           ec: ExecutionContext) = {
    this(DefaultStandardExpressionFunctionsParams(pathBuilders, callContext, ioActorEndpoint, ec))
  }

  override def postMapping(path: Path) = {
    path match {
      case _: DefaultPath if !path.isAbsolute => callContext.root.resolve(path)
      case _ => path
    }
  }
}
