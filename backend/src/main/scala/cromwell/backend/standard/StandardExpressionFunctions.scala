package cromwell.backend.standard

import akka.actor.ActorRef
import cromwell.backend.io.GlobFunctions
import cromwell.backend.wdl.{ReadLikeFunctions, WriteFunctions}
import cromwell.core.CallContext
import cromwell.core.io.{AsyncIo, DefaultIoCommandBuilder, IoCommandBuilder}
import cromwell.core.path.{Path, PathBuilder}
import wom.values.{WomSingleFile, WomValue}

import scala.concurrent.ExecutionContext
import scala.util.{Success, Try}

trait StandardExpressionFunctionsParams {
  def pathBuilders: List[PathBuilder]

  def callContext: CallContext

  def ioActorEndpoint: ActorRef
  
  def executionContext: ExecutionContext
}

case class DefaultStandardExpressionFunctionsParams(override val pathBuilders: List[PathBuilder],
                                                    override val callContext: CallContext,
                                                    override val ioActorEndpoint: ActorRef,
                                                    override val executionContext: ExecutionContext
                                                   ) extends StandardExpressionFunctionsParams

// TODO: Once we figure out premapping and postmapping, maybe we can standardize that behavior. Currently that's the most important feature that subclasses override.
class StandardExpressionFunctions(val standardParams: StandardExpressionFunctionsParams)
  extends GlobFunctions with ReadLikeFunctions with WriteFunctions {

  override lazy val ec = standardParams.executionContext
  
  protected lazy val ioCommandBuilder: IoCommandBuilder = DefaultIoCommandBuilder

  override lazy val asyncIo = new AsyncIo(standardParams.ioActorEndpoint, ioCommandBuilder)

  override val pathBuilders: List[PathBuilder] = standardParams.pathBuilders

  val callContext: CallContext = standardParams.callContext

  val writeDirectory: Path = callContext.root

  override def stdout(params: Seq[Try[WomValue]]): Try[WomSingleFile] = Success(WomSingleFile(callContext.stdout))

  override def stderr(params: Seq[Try[WomValue]]): Try[WomSingleFile] = Success(WomSingleFile(callContext.stderr))
}
