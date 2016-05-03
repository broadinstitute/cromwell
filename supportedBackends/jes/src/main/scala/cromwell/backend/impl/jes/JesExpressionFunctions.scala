package cromwell.backend.impl.jes

import java.nio.file.{FileSystem, Path}

import better.files._
import cromwell.backend.wdl.{PureFunctions, ReadLikeFunctions, WriteFunctions}
import cromwell.core.CallContext
import wdl4s.expression.WdlStandardLibraryFunctions
import wdl4s.values._

import scala.language.postfixOps
import scala.util.{Success, Try}

class JesExpressionFunctions(override val fileSystems: List[FileSystem],
                             context: CallContext
                                ) extends WdlStandardLibraryFunctions with PureFunctions with ReadLikeFunctions with WriteFunctions {

  private def globDirectory(glob: String): String = s"glob-${glob.md5Sum}/"

  override def globPath(glob: String): String = context.root.resolve(globDirectory(glob)).toString

  override def glob(path: String, pattern: String): Seq[String] = {
    path.toAbsolutePath(fileSystems).asDirectory.glob("**/*") map { _.path.fullPath } filterNot { _.toString == path } toSeq
  }

  override def postMapping(path: Path): Path = if (!path.isAbsolute) context.root.resolve(path) else path

  override def stdout(params: Seq[Try[WdlValue]]) = Success(WdlFile(context.stdout))
  override def stderr(params: Seq[Try[WdlValue]]) = Success(WdlFile(context.stderr))

  override val writeDirectory: Path = context.root
}
