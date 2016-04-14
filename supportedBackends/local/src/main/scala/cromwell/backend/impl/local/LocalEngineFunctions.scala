package cromwell.backend.impl.local

import java.nio.file.{FileSystem, Path}

import cromwell.core.{CallContext, WorkflowContext}
import cromwell.{CallEngineFunctions, WorkflowEngineFunctions}
import wdl4s.values.WdlValue

import scala.language.postfixOps
import scala.util.Try

object LocalCallEngineFunctions {
  private val LocalFSScheme = "file"
  def isLocalPath(path: Path) = path.toUri.getScheme == LocalCallEngineFunctions.LocalFSScheme
}

class LocalWorkflowEngineFunctions(fileSystems: List[FileSystem], context: WorkflowContext) extends WorkflowEngineFunctions(context) {
  import better.files._
  import cromwell.core.PathFactory._
  import LocalCallEngineFunctions._

  override def globPath(glob: String): String = context.root
  override def glob(path: String, pattern: String): Seq[String] = {
    toPath(path).glob(s"**/$pattern") map { _.path.fullPath } toSeq
  }

  override def toPath(str: String): Path = {
    val path = buildPath(str, fileSystems)
    if (!path.isAbsolute && isLocalPath(path)) buildPath(context.root, fileSystems).resolve(path) else path
  }
}

class LocalCallEngineFunctions(fileSystems: List[FileSystem], context: CallContext) extends LocalWorkflowEngineFunctions(fileSystems ,context) with CallEngineFunctions {
  override def stdout(params: Seq[Try[WdlValue]]) = stdout(context)
  override def stderr(params: Seq[Try[WdlValue]]) = stderr(context)
}

