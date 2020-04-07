package cromwell.backend.sfs

import akka.actor.ActorRef
import cromwell.backend.io._
import cromwell.backend.standard.{DefaultStandardExpressionFunctionsParams, StandardExpressionFunctions, StandardExpressionFunctionsParams}
import cromwell.core.CallContext
import cromwell.core.path.{DefaultPath, DefaultPathBuilder, Path, PathBuilder, PathFactory}

import scala.concurrent.ExecutionContext

object SharedFileSystemExpressionFunctions {
  def apply(jobPaths: JobPaths,
            pathBuilders: List[PathBuilder],
            ioActorProxy: ActorRef,
            ec: ExecutionContext,
            forInput: Boolean): SharedFileSystemExpressionFunctions = {
    new SharedFileSystemExpressionFunctions(pathBuilders, jobPaths.callContext, ioActorProxy, ec, forInput)
  }
}

class SharedFileSystemExpressionFunctions(standardParams: StandardExpressionFunctionsParams)
  extends StandardExpressionFunctions(standardParams) {

  def this(pathBuilders: List[PathBuilder],
           callContext: CallContext,
           ioActorProxy: ActorRef,
           ec: ExecutionContext,
           forInput: Boolean) = {
    this(DefaultStandardExpressionFunctionsParams(pathBuilders, callContext, ioActorProxy, ec, forInput))
  }

  lazy val cromwellCwd: Path = DefaultPathBuilder.build(new java.io.File(".").getCanonicalPath).get
  def noCallContextPostMapping(path: Path): Path = {
    path match {
      case _: DefaultPath if !path.isAbsolute => cromwellCwd.resolve(path)
      case _ => path
    }
  }

  override def postMapping(path: Path) = {
    path match {
      case _: DefaultPath if !path.isAbsolute => callContext.root.resolve(path)
      case _ => path
    }
  }

  override def buildPath(string: String): Path = PathFactory.buildPath(string, pathBuilders, preMapping, if (forInput) noCallContextPostMapping else postMapping)
}
