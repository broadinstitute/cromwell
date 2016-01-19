package cromwell.backend.provider.local

import java.nio.file.{Paths, Path}

import cromwell.backend.expressions.WdlStandardLibraryImpl
import cromwell.backend.io.IoInterface
import cromwell.backend.io.shared.SharedFileSystemIoInterface
import cromwell.backend.provider.local.LocalBackend._
import wdl4s.values.WdlFile

import scala.util.{Success, Try}

class WorkflowEngineFunctions(val executionDir: Path) extends WdlStandardLibraryImpl {
  override val interface: IoInterface = SharedFileSystemIoInterface.instance
  override def tempFilePath: String = executionDir.toString
  override def globPath(glob: String): String = executionDir.toString
}

trait CallEngineFunctions extends WdlStandardLibraryImpl {
  protected def stdout(executionDir: Path): Try[WdlFile] = Success(WdlFile(Paths.get(executionDir.toString, StdoutFile).toString))
  protected def stderr(executionDir: Path): Try[WdlFile] = Success(WdlFile(Paths.get(executionDir.toString, StderrFile).toString))
}