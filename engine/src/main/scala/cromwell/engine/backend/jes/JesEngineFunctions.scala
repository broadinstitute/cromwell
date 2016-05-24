package cromwell.engine.backend.jes

import java.nio.file.FileSystem

import better.files._
import cromwell.engine.backend.io._
import cromwell.engine.backend.{CallContext, CallEngineFunctions, WorkflowContext, WorkflowEngineFunctions}
import wdl4s.values._

import scala.language.postfixOps
import scala.util.Try
class JesWorkflowEngineFunctions(fileSystems: List[FileSystem], context: WorkflowContext) extends WorkflowEngineFunctions(fileSystems, context) {
  override def globPath(glob: String): String = context.root.toAbsolutePath(fileSystems).resolve(JesBackend.globDirectory(glob)).toString
  override def glob(path: String, pattern: String): Seq[String] = {
    path.toAbsolutePath(fileSystems).asDirectory.list map { _.path.fullPath } toSeq
  }
  override def adjustFilePath(path: String) = if (!path.isGcsUrl) context.root.toAbsolutePath(fileSystems).resolve(path).toString else path
}

class JesCallEngineFunctions(fileSystems: List[FileSystem], context: CallContext) extends JesWorkflowEngineFunctions(fileSystems, context) with CallEngineFunctions {
  override def stdout(params: Seq[Try[WdlValue]]) = stdout(context)
  override def stderr(params: Seq[Try[WdlValue]]) = stderr(context)
}
