package cromwell.backend.impl.jes

import java.nio.file.{FileSystem, Path}

import better.files._
import cromwell.backend.wdl.{PureFunctions, ReadLikeFunctions, WriteFunctions}
import cromwell.backend.impl.jes.JesImplicits.PathString
import cromwell.core.CallContext
import cromwell.filesystems.gcs.GcsFileSystem
import wdl4s.expression.WdlStandardLibraryFunctions
import wdl4s.values._

import scala.language.postfixOps
import scala.util.{Success, Try}

class JesExpressionFunctions(override val fileSystems: List[FileSystem],
                             context: CallContext
                             ) extends WdlStandardLibraryFunctions with PureFunctions with ReadLikeFunctions with WriteFunctions {
  import JesExpressionFunctions.EnhancedPath

  private def globDirectory(glob: String): String = s"glob-${glob.md5Sum}/"

  override def globPath(glob: String): String = context.root.resolve(globDirectory(glob)).toString

  override def glob(path: String, pattern: String): Seq[String] = {
    path.toAbsolutePath(fileSystems).asDirectory.glob("**/*") map { _.path.fullPath } filterNot { _.toString == path } toSeq
  }

  override def preMapping(str: String): String = if (!GcsFileSystem.isAbsoluteGcsPath(str)) context.root.resolve(str).toString else str

  override def stdout(params: Seq[Try[WdlValue]]) = Success(WdlFile(context.stdout))
  override def stderr(params: Seq[Try[WdlValue]]) = Success(WdlFile(context.stderr))

  override val writeDirectory: Path = context.root
}

object JesExpressionFunctions {
  implicit class EnhancedPath(val path: Path) extends AnyVal {
    def asDirectory = path.toString.toDirectory(path.getFileSystem)
  }
}
