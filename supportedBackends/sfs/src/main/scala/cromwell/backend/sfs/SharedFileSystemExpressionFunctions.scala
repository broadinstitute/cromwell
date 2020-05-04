package cromwell.backend.sfs

import akka.actor.ActorRef
import cromwell.backend.io._
import cromwell.backend.standard.{DefaultStandardExpressionFunctionsParams, StandardExpressionFunctions, StandardExpressionFunctionsParams}
import cromwell.core.CallContext
import cromwell.core.path.{DefaultPath, DefaultPathBuilder, Path, PathBuilder}
import wom.expression.IoFunctionSet

import scala.concurrent.ExecutionContext

object SharedFileSystemExpressionFunctions {
  def apply(jobPaths: JobPaths,
            pathBuilders: List[PathBuilder],
            ioActorProxy: ActorRef,
            ec: ExecutionContext): SharedFileSystemExpressionFunctions = {
    new SharedFileSystemExpressionFunctions(pathBuilders, jobPaths.callContext, ioActorProxy, ec)
  }
}

class SharedFileSystemExpressionFunctions(standardParams: StandardExpressionFunctionsParams)
  extends StandardExpressionFunctions(standardParams) {

  def this(pathBuilders: List[PathBuilder],
           callContext: CallContext,
           ioActorProxy: ActorRef,
           ec: ExecutionContext) = {
    this(DefaultStandardExpressionFunctionsParams(pathBuilders, callContext, ioActorProxy, ec))
  }

  // Use a separate class for forInput = true, otherwise makeInputSpecificFunctions is self-referential and will fail.
  override def makeInputSpecificFunctions: IoFunctionSet = new SharedFileSystemExpressionFunctionsForInput(standardParams)

  lazy val forInput = false
  lazy val cromwellCwd: Path = DefaultPathBuilder.build(sys.props("user.dir")).get

  override def postMapping(path: Path) = {
    path match {
      case _: DefaultPath if !path.isAbsolute && forInput => cromwellCwd.resolve(path)
      case _: DefaultPath if !path.isAbsolute => callContext.root.resolve(path)
      case _ => path
    }
  }
}

class SharedFileSystemExpressionFunctionsForInput(standardParams: StandardExpressionFunctionsParams)
  extends SharedFileSystemExpressionFunctions(standardParams) {
  override def makeInputSpecificFunctions: IoFunctionSet = this
  override lazy val forInput = true
}
