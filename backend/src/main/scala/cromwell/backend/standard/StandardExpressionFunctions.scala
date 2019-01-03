package cromwell.backend.standard

import akka.actor.ActorRef
import cromwell.backend.io.{DirectoryFunctions, GlobFunctions}
import cromwell.backend.{ReadLikeFunctions, WriteFunctions}
import cromwell.core.CallContext
import cromwell.core.io._
import cromwell.core.path.PathFactory.PathBuilders
import cromwell.core.path.{Path, PathBuilder}

import scala.concurrent.ExecutionContext

trait StandardExpressionFunctionsParams {
  def pathBuilders: PathBuilders

  def callContext: CallContext

  def ioActorProxy: ActorRef

  def executionContext: ExecutionContext
}

case class DefaultStandardExpressionFunctionsParams(override val pathBuilders: PathBuilders,
                                                    override val callContext: CallContext,
                                                    override val ioActorProxy: ActorRef,
                                                    override val executionContext: ExecutionContext
                                                   ) extends StandardExpressionFunctionsParams

// TODO: Once we figure out premapping and postmapping, maybe we can standardize that behavior. Currently that's the most important feature that subclasses override.
class StandardExpressionFunctions(val standardParams: StandardExpressionFunctionsParams)
  extends GlobFunctions with DirectoryFunctions with ReadLikeFunctions with WriteFunctions with CallCorePathFunctions {

  override lazy val ec = standardParams.executionContext

  protected lazy val ioCommandBuilder: IoCommandBuilder = DefaultIoCommandBuilder

  override lazy val asyncIo = new AsyncIo(standardParams.ioActorProxy, ioCommandBuilder)

  override val pathBuilders: List[PathBuilder] = standardParams.pathBuilders

  val callContext: CallContext = standardParams.callContext

  val writeDirectory: Path = callContext.root
  
  val isDocker: Boolean = callContext.isDocker
}
