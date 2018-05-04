package cromwell.engine

import cromwell.backend.wdl.ReadLikeFunctions
import cromwell.core.io.{AsyncIo, WorkflowCorePathFunctions}
import cromwell.core.path.PathBuilder
import wom.values.WomSingleFile

import scala.concurrent.{ExecutionContext, Future}

class EngineIoFunctions(val pathBuilders: List[PathBuilder], override val asyncIo: AsyncIo, override val ec: ExecutionContext) extends ReadLikeFunctions with WorkflowCorePathFunctions {
  override def glob(pattern: String): Future[Seq[String]] = throw new NotImplementedError(s"glob(path, pattern) not implemented yet")

  override def writeFile(path: String, content: String): Future[WomSingleFile] =
    Future.failed(new Exception("Cromwell does not support writing files from a workflow context"))

  override def copyFile(pathFrom: String, targetName: String): Future[WomSingleFile] =
    Future.failed(new Exception("Cromwell does not support copying files from a workflow context"))

  override def listDirectory(path: String)(visited: Vector[String]) = throw new NotImplementedError(s"listDirectory not implemented yet")

  override def isDirectory(path: String) = Future.successful(buildPath(path).isDirectory)
}
