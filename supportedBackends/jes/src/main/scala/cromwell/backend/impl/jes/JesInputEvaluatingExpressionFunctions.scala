package cromwell.backend.impl.jes

import java.nio.file.Path

import cromwell.backend.impl.jes.JesAsyncBackendJobExecutionActor.GlobComplete
import cromwell.backend.wdl.{PureFunctions, ReadLikeFunctions, WriteFunctions}
import cromwell.core.CallContext
import cromwell.core.path.PathBuilder
import cromwell.filesystems.gcs.GcsPathBuilder
import wdl4s.expression.WdlStandardLibraryFunctions
import wdl4s.values._

import scala.util.{Success, Try}

class JesInputEvaluatingExpressionFunctions(override val pathBuilders: List[PathBuilder], val context: CallContext)
  extends WdlStandardLibraryFunctions with PureFunctions with ReadLikeFunctions with WriteFunctions {

  private def globDirectory(glob: String): String = s"glob-${glob.md5Sum}/"
  override def globPath(glob: String): String = context.root.resolve(globDirectory(glob)).toString

  override def glob(path: String, pattern: String): Seq[String] = {
    throw new RuntimeException("The WDL function glob() should only be used in a task's output block")
  }

  override def preMapping(str: String): String = if (!GcsPathBuilder.isValidGcsUrl(str)) {
    context.root.resolve(str.stripPrefix("/")).toUri.toString
  } else str

  override def stdout(params: Seq[Try[WdlValue]]) = Success(WdlFile(context.stdout))
  override def stderr(params: Seq[Try[WdlValue]]) = Success(WdlFile(context.stderr))

  override val writeDirectory: Path = context.root
}

class JesOutputEvaluatingExpressionFunctions(pathBuilders: List[PathBuilder], context: CallContext, evaluatedGlobs: Seq[GlobComplete])
  extends JesInputEvaluatingExpressionFunctions(pathBuilders, context) {

  override def glob(path: String, pattern: String): Seq[String] = {

    val directory = context.root.resolve(s"glob-${pattern.md5Sum}/").toRealPath().toUri

    evaluatedGlobs find { globComplete: GlobComplete => globComplete.globInfo.path.equals(directory.toString) } match {
      case Some(globComplete) => globComplete.globContents
      case None => throw new RuntimeException(s"Path '$path' was not declared by JES as a glob output of the call")
    }
  }
}
