package cromwell.backend.impl.jes

import java.nio.file.{Files, Path}

import cromwell.backend.wdl.{PureFunctions, ReadLikeFunctions, WriteFunctions}
import cromwell.core.CallContext
import cromwell.core.path.PathBuilder
import cromwell.filesystems.gcs.GcsPathBuilder
import wdl4s.expression.WdlStandardLibraryFunctions
import wdl4s.values._

import scala.collection.JavaConverters._
import scala.language.postfixOps
import scala.util.{Success, Try}

class JesExpressionFunctions(override val pathBuilders: List[PathBuilder], context: CallContext)
  extends WdlStandardLibraryFunctions with PureFunctions with ReadLikeFunctions with WriteFunctions {

  private def globDirectory(glob: String): String = s"glob-${glob.md5Sum}/"

  override def globPath(glob: String): String = context.root.resolve(globDirectory(glob)).toString

  override def glob(path: String, pattern: String): Seq[String] = {
    /* Globbing is not implemented in the gcs-nio implementation yet (specifically PathMatcher) at this time (09/2016)
     * which is fine here since all files in the globing directory have already be matched against the glob pattern by JES.
     *
     * Also, better.file.File can't be used here because it calls toAbsolutePath on the path in File.apply, which adds a leading "/" in the path.
     * This makes the newDirectoryStream implementation of CloudStorageFileSystemProvider return nothing
     * because it doesn't call toRealPath for this method before doing I/O and hence does not remove the "/" prefix (bug ?)
     * See com.google.cloud.storage.contrib.nio.CloudStorageConfiguration#stripPrefixSlash()
     */
    val directory = context.root.resolve(s"glob-${pattern.md5Sum}/").toRealPath()
    Files.list(directory).iterator().asScala filterNot { Files.isDirectory(_) } map { _.toUri.toString } toSeq
  }

  override def preMapping(str: String): String = if (!GcsPathBuilder.isValidGcsUrl(str)) {
    context.root.resolve(str.stripPrefix("/")).toUri.toString
  } else str

  override def stdout(params: Seq[Try[WdlValue]]) = Success(WdlFile(context.stdout))
  override def stderr(params: Seq[Try[WdlValue]]) = Success(WdlFile(context.stderr))

  override val writeDirectory: Path = context.root
}
