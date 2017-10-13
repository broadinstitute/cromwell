package cromwell.backend.standard

import cromwell.backend.io.GlobFunctions
import cromwell.backend.wdl.{ReadLikeFunctions, WriteFunctions}
import cromwell.core.CallContext
import cromwell.core.path.{Path, PathBuilder}
import wom.values.{WdlFile, WdlValue}

import scala.util.{Success, Try}

trait StandardExpressionFunctionsParams {
  def pathBuilders: List[PathBuilder]

  def callContext: CallContext
}

case class DefaultStandardExpressionFunctionsParams(override val pathBuilders: List[PathBuilder],
                                                    override val callContext: CallContext
                                                   ) extends StandardExpressionFunctionsParams

// TODO: Once we figure out premapping and postmapping, maybe we can standardize that behavior. Currently that's the most important feature that subclasses override.
class StandardExpressionFunctions(val standardParams: StandardExpressionFunctionsParams)
  extends GlobFunctions with ReadLikeFunctions with WriteFunctions {

  override val pathBuilders: List[PathBuilder] = standardParams.pathBuilders

  val callContext: CallContext = standardParams.callContext

  val writeDirectory: Path = callContext.root

  override def stdout(params: Seq[Try[WdlValue]]): Try[WdlFile] = Success(WdlFile(callContext.stdout))

  override def stderr(params: Seq[Try[WdlValue]]): Try[WdlFile] = Success(WdlFile(callContext.stderr))
}
