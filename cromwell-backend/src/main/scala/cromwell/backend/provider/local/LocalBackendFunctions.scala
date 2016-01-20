package cromwell.backend.provider.local

import java.nio.file.{Path, Paths}

import cromwell.backend.expressions.WdlStandardLibraryImpl
import cromwell.backend.io.IoInterface
import cromwell.backend.io.shared.SharedFileSystemIoInterface
import cromwell.backend.provider.local.LocalBackend._
import wdl4s.values.WdlFile

import scala.util.{Success, Try}

class WorkflowEngineFunctions(val executionDir: Path) extends WdlStandardLibraryImpl {

  def isUriWithProtocol(str: String): Boolean = "^[a-z]+://".r.findFirstIn(str).nonEmpty

  override val interface: IoInterface = SharedFileSystemIoInterface.instance

  override def tempFilePath: String = executionDir.toString

  override def globPath(glob: String): String = executionDir.toString

  protected def stdout(executionDir: Path): Try[WdlFile] = Success(WdlFile(Paths.get(executionDir.toString, StdoutFile).toString))

  protected def stderr(executionDir: Path): Try[WdlFile] = Success(WdlFile(Paths.get(executionDir.toString, StderrFile).toString))

  override def fileContentsToString(path: String) = {
    if (!Paths.get(path).isAbsolute && !isUriWithProtocol(path))
      interface.readFile(executionDir.resolve(path).toString)
    else interface.readFile(path)
  }
}