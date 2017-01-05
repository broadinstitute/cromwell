package cromwell.backend.impl.jes

import java.nio.file.Path

import cromwell.backend.io.GlobFunctions
import cromwell.backend.wdl.{ReadLikeFunctions, WriteFunctions}
import cromwell.core.CallContext
import cromwell.core.path.PathBuilder
import cromwell.filesystems.gcs.GcsPathBuilder
import cromwell.core.path.PathImplicits._
import wdl4s.expression.{PureStandardLibraryFunctionsLike, WdlStandardLibraryFunctions}
import wdl4s.values._

import scala.util.{Success, Try}

class JesExpressionFunctions(override val pathBuilders: List[PathBuilder], context: CallContext)
  extends WdlStandardLibraryFunctions with PureStandardLibraryFunctionsLike with ReadLikeFunctions with WriteFunctions with GlobFunctions {

  def callContext: CallContext = context

  override def writeTempFile(path: String, prefix: String, suffix: String, content: String): String = super[WriteFunctions].writeTempFile(path, prefix, suffix, content)

  override def preMapping(str: String): String = if (!GcsPathBuilder.isValidGcsUrl(str)) {
    context.root.resolve(str.stripPrefix("/")).toRealString
  } else str

  override def stdout(params: Seq[Try[WdlValue]]) = Success(WdlFile(context.stdout))
  override def stderr(params: Seq[Try[WdlValue]]) = Success(WdlFile(context.stderr))

  override val writeDirectory: Path = context.root
}
