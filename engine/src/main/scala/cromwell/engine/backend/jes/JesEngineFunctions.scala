package cromwell.engine.backend.jes

import java.nio.file.{FileSystem, Path}

import better.files._
import cromwell.backend.wdl.{OldCallEngineFunctions, OldWorkflowEngineFunctions}
import cromwell.core.{OldCallContext, OldWorkflowContext, PathFactory}
import cromwell.engine.backend.io._
import wdl4s.values._

import scala.language.postfixOps
import scala.util.Try
@deprecated(message = "This class will not be part of the PBE universe", since = "May 2nd 2016")
class OldJesWorkflowEngineFunctions(val fileSystems: List[FileSystem], context: OldWorkflowContext) extends OldWorkflowEngineFunctions(context) with PathFactory {
  override def globPath(glob: String): String = context.root.toAbsolutePath(fileSystems).resolve(OldStyleJesBackend.globDirectory(glob)).toString
  override def glob(path: String, pattern: String): Seq[String] = {
    path.toAbsolutePath(fileSystems).asDirectory.glob("**/*") map { _.path.fullPath } filterNot { _.toString == path } toSeq
  }

  override def postMapping(path: Path) = if (!path.isAbsolute) context.root.toAbsolutePath(fileSystems).resolve(path) else path
}

@deprecated(message = "This class will not be part of the PBE universe", since = "May 2nd 2016")
class OldJesCallEngineFunctions(override val fileSystems: List[FileSystem], context: OldCallContext) extends OldJesWorkflowEngineFunctions(fileSystems, context) with OldCallEngineFunctions {
  override def stdout(params: Seq[Try[WdlValue]]) = stdout(context)
  override def stderr(params: Seq[Try[WdlValue]]) = stderr(context)
}
