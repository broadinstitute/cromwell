package cromwell.core

import wom.expression.IoFunctionSet
import wom.values.{WomSingleFile, WomValue}

import scala.concurrent.Future
import scala.util.{Failure, Try}

case object NoIoFunctionSet extends IoFunctionSet {
  override def readFile(path: String, maxBytes: Option[Int] = None, failOnOverflow: Boolean = false): Future[String] = Future.failed(new NotImplementedError("readFile is not available here"))

  override def writeFile(path: String, content: String): Future[WomSingleFile] = {
    Future.failed(new NotImplementedError("writeFile is not available here"))
  }

  override def copyFile(pathFrom: String, targetName: String): Future[WomSingleFile] = {
    throw new Exception("copyFile is not available here")
  }

  override def stdout(params: Seq[Try[WomValue]]): Try[WomSingleFile] = {
    Failure(new NotImplementedError("stdout is not available here"))
  }

  override def stderr(params: Seq[Try[WomValue]]): Try[WomSingleFile] = {
    Failure(new NotImplementedError("stderr is not available here"))
  }

  override def glob(pattern: String): Future[Seq[String]] = throw new NotImplementedError("glob is not available here")

  override def listAllFilesUnderDirectory(dirPath: String): Nothing =
    throw new NotImplementedError("listAllFilesUnderDirectory is not available here")

  override def size(file: WomValue): Future[Double] = Future.failed(new NotImplementedError("size is not available here"))
}
