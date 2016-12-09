package cromwell.backend.impl.jes

import java.nio.file.{Files, Path}

import cromwell.backend.wdl.{ReadLikeFunctions, WriteFunctions}
import cromwell.core.CallContext
import cromwell.core.path.PathBuilder
import cromwell.filesystems.gcs.GcsPathBuilder
import cromwell.core.path.PathImplicits._
import wdl4s.expression.{PureStandardLibraryFunctionsLike, WdlStandardLibraryFunctions}
import wdl4s.values._

import scala.collection.JavaConverters._
import scala.util.{Success, Try}

class JesExpressionFunctions(override val pathBuilders: List[PathBuilder], context: CallContext)
  extends WdlStandardLibraryFunctions with PureStandardLibraryFunctionsLike with ReadLikeFunctions with WriteFunctions {

  override def writeTempFile(path: String, prefix: String, suffix: String, content: String): String = super[WriteFunctions].writeTempFile(path, prefix, suffix, content)
  private[jes] def globDirectory(glob: String): String = globName(glob) + "/"
  private[jes] def globName(glob: String) = s"glob-${glob.md5Sum}"

  override def globPath(glob: String): String = context.root.resolve(globDirectory(glob)).toString

  override def glob(path: String, pattern: String): Seq[String] = {
    val name = globName(pattern)
    val listFile = context.root.resolve(s"$name.list").toRealPath()
    Files.readAllLines(listFile).asScala map { fileName => context.root.resolve(s"$name/$fileName").toRealString }
  }

  override def preMapping(str: String): String = if (!GcsPathBuilder.isValidGcsUrl(str)) {
    context.root.resolve(str.stripPrefix("/")).toRealString
  } else str

  override def stdout(params: Seq[Try[WdlValue]]) = Success(WdlFile(context.stdout))
  override def stderr(params: Seq[Try[WdlValue]]) = Success(WdlFile(context.stderr))

  override val writeDirectory: Path = context.root
}
