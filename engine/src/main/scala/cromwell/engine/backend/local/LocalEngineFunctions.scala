package cromwell.engine.backend.local

import java.nio.file.{FileSystem, Path}

import cromwell.backend.wdl.{OldCallEngineFunctions, OldWorkflowEngineFunctions}
import cromwell.core.{OldCallContext, OldWorkflowContext, PathFactory}
import cromwell.engine.backend
import wdl4s.values.WdlValue

import scala.language.postfixOps
import scala.util.Try

@deprecated(message = "This class will not be part of the PBE universe", since = "May 2nd 2016")
object OldStyleLocalWorkflowEngineFunctions {
  private val LocalFSScheme = "file"
  def isLocalPath(path: Path) = path.toUri.getScheme == LocalFSScheme
}

@deprecated(message = "This class will not be part of the PBE universe", since = "May 2nd 2016")
class OldStyleLocalWorkflowEngineFunctions(val fileSystems: List[FileSystem], context: OldWorkflowContext) extends OldWorkflowEngineFunctions(context) with PathFactory {
  import OldStyleLocalWorkflowEngineFunctions._
  import backend.io._
  import better.files._

  override def globPath(glob: String): String = context.root
  override def glob(path: String, pattern: String): Seq[String] = {
    path.toAbsolutePath(fileSystems).glob(s"**/$pattern") map { _.path.fullPath } toSeq
  }

  override def postMapping(path: Path) = if (!path.isAbsolute && isLocalPath(path)) context.root.toPath(fileSystems).resolve(path) else path
}

@deprecated(message = "This class will not be part of the PBE universe", since = "May 2nd 2016")
class OldStyleLocalCallEngineFunctions(fileSystems: List[FileSystem], context: OldCallContext) extends OldStyleLocalWorkflowEngineFunctions(fileSystems ,context) with OldCallEngineFunctions {
  override def stdout(params: Seq[Try[WdlValue]]) = stdout(context)
  override def stderr(params: Seq[Try[WdlValue]]) = stderr(context)
}

