package cromwell.engine

import cromwell.backend.ReadLikeFunctions
import cromwell.core.io.{AsyncIo, WorkflowCorePathFunctions}
import cromwell.core.path.PathBuilder
import wom.values.WomSingleFile
import better.files.File._

import scala.concurrent.{ExecutionContext, Future}

class EngineIoFunctions(val pathBuilders: List[PathBuilder], override val asyncIo: AsyncIo, override val ec: ExecutionContext) extends ReadLikeFunctions with WorkflowCorePathFunctions {
  override def glob(pattern: String): Future[Seq[String]] = throw new UnsupportedOperationException(s"glob(path, pattern) not implemented yet")

  // TODO: This is not suited for multi backend / multi filesystem use. Keep local for now to not break local CWL conf tests
  override def writeFile(path: String, content: String): Future[WomSingleFile] = Future.successful {
    val cromwellPath = buildPath(path)
    val string = if (cromwellPath.isAbsolute) 
      cromwellPath.write(content).pathAsString
    else 
      (newTemporaryDirectory() / path).write(content).pathAsString
    WomSingleFile(string)
  }

  override def copyFile(pathFrom: String, targetName: String): Future[WomSingleFile] =
    Future.failed(new Exception("Cromwell does not support copying files from a workflow context"))

  override def listAllFilesUnderDirectory(dirPath: String): Nothing =
    throw new UnsupportedOperationException(s"listAllFilesUnderDirectory not implemented yet")

  override def listDirectory(path: String)(visited: Vector[String]) = throw new UnsupportedOperationException(s"listDirectory not implemented yet")

  override def isDirectory(path: String) = Future.successful(buildPath(path).isDirectory)

  // TODO: This is not suited for multi backend / multi filesystem use. Keep local for now to not break local CWL conf tests
  override def createTemporaryDirectory(name: Option[String]) = Future.successful {
    name map {
      newTemporaryDirectory().createChild(_, asDirectory = true).pathAsString
    } getOrElse newTemporaryDirectory().pathAsString
  }
}
