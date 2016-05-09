package cromwell.engine.backend.local

import java.nio.file.{FileSystem, Path}

import cromwell.core.{CallContext, WorkflowContext}
import cromwell.engine.backend
import cromwell.{core, CallEngineFunctions, WorkflowEngineFunctions}
import wdl4s.values.WdlValue

import scala.language.postfixOps
import scala.util.Try

@deprecated(message = "This class will not be part of the PBE universe", since = "May 2nd 2016")
object OldStyleLocalWorkflowEngineFunctions {
  private val LocalFSScheme = "file"
  def isLocalPath(path: Path) = path.toUri.getScheme == LocalFSScheme
}

@deprecated(message = "This class will not be part of the PBE universe", since = "May 2nd 2016")
class OldStyleLocalWorkflowEngineFunctions(fileSystems: List[FileSystem], context: WorkflowContext) extends WorkflowEngineFunctions(context) {
  import backend.io._
  import core.PathFactory._
  import better.files._
  import OldStyleLocalWorkflowEngineFunctions._

  override def globPath(glob: String): String = context.root
  override def glob(path: String, pattern: String): Seq[String] = {
    path.toAbsolutePath(fileSystems).glob(s"**/$pattern") map { _.path.fullPath } toSeq
  }

  override def toPath(str: String): Path = {
    val path = buildPath(str, fileSystems)
    if (!path.isAbsolute && isLocalPath(path)) context.root.toPath(fileSystems).resolve(path) else path
  }
}

@deprecated(message = "This class will not be part of the PBE universe", since = "May 2nd 2016")
class OldStyleLocalCallEngineFunctions(fileSystems: List[FileSystem], context: CallContext) extends OldStyleLocalWorkflowEngineFunctions(fileSystems ,context) with CallEngineFunctions {
  override def stdout(params: Seq[Try[WdlValue]]) = stdout(context)
  override def stderr(params: Seq[Try[WdlValue]]) = stderr(context)
}

