package cromwell.engine

import cromwell.backend.wdl.ReadLikeFunctions
import cromwell.core.io.AsyncIo
import cromwell.core.path.PathBuilder
import wom.values.{WomSingleFile, WomValue}

import scala.concurrent.{ExecutionContext, Future}
import scala.util.{Failure, Try}

class EngineIoFunctions(val pathBuilders: List[PathBuilder], override val asyncIo: AsyncIo, override val ec: ExecutionContext) extends ReadLikeFunctions {
  private def fail(name: String) = Failure(new NotImplementedError(s"$name() not supported at the workflow level yet"))

  override def stdout(params: Seq[Try[WomValue]]): Try[WomSingleFile] = fail("stdout")
  override def stderr(params: Seq[Try[WomValue]]): Try[WomSingleFile] = fail("stderr")
  override def glob(pattern: String): Future[Seq[String]] = throw new NotImplementedError(s"glob(path, pattern) not implemented yet")

  override def writeFile(path: String, content: String): Future[WomSingleFile] =
    Future.failed(new Exception("Cromwell does not support writing files from a workflow context"))

  override def copyFile(pathFrom: String, targetName: String): Future[WomSingleFile] =
    Future.failed(new Exception("Cromwell does not support copying files from a workflow context"))

  override def listAllFilesUnderDirectory(dirPath: String): Nothing =
    throw new NotImplementedError(s"listAllFilesUnderDirectory not implemented yet")
}
