package cromwell.engine.backend.jes

import java.nio.file.{FileSystem, Path}

import better.files._
import cromwell.core.{CallContext, WorkflowContext}
import cromwell.engine.backend.io._
import cromwell.{CallEngineFunctions, WorkflowEngineFunctions}
import wdl4s.values._

import scala.language.postfixOps
import scala.util.Try
class JesWorkflowEngineFunctions(fileSystems: List[FileSystem], context: WorkflowContext) extends WorkflowEngineFunctions(context) {
  import cromwell.core.PathFactory._

  override def globPath(glob: String): String = context.root.toAbsolutePath(fileSystems).resolve(JesBackend.globDirectory(glob)).toString
  override def glob(path: String, pattern: String): Seq[String] = {
    path.toAbsolutePath(fileSystems).asDirectory.glob("**/*") map { _.path.fullPath } filterNot { _.toString == path } toSeq
  }

  override def toPath(str: String): Path = {
    val path = buildPath(str, fileSystems)
    if (!path.isAbsolute) context.root.toAbsolutePath(fileSystems).resolve(path) else path
  }
}

class JesCallEngineFunctions(fileSystems: List[FileSystem], context: CallContext) extends JesWorkflowEngineFunctions(fileSystems, context) with CallEngineFunctions {
  override def stdout(params: Seq[Try[WdlValue]]) = stdout(context)
  override def stderr(params: Seq[Try[WdlValue]]) = stderr(context)
}
