package cromwell.engine.backend.local

import java.nio.file.FileSystem

import cromwell.engine.{CallContext, CallEngineFunctions, WorkflowContext, WorkflowEngineFunctions, _}
import wdl4s.values.WdlValue

import scala.language.postfixOps
import scala.util.Try

class LocalWorkflowEngineFunctions(fileSystems: List[FileSystem], context: WorkflowContext) extends WorkflowEngineFunctions(fileSystems, context) {
  import backend.io._
  import better.files._

  override def globPath(glob: String): String = context.root
  override def glob(path: String, pattern: String): Seq[String] = {
    path.toAbsolutePath(fileSystems).glob(s"**/$pattern") map { _.path.fullPath } toSeq
  }
}

class LocalCallEngineFunctions(fileSystems: List[FileSystem], context: CallContext) extends LocalWorkflowEngineFunctions(fileSystems ,context) with CallEngineFunctions {
  import backend.io._

  override def adjustFilePath(path: String) = {
    if (!path.toPath(fileSystems).isAbsolute && !path.isUriWithProtocol) context.root.toPath(fileSystems).resolve(path).toString else path
  }

  override def stdout(params: Seq[Try[WdlValue]]) = stdout(context)
  override def stderr(params: Seq[Try[WdlValue]]) = stderr(context)
}

