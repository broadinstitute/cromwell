package cromwell.core

import wom.expression.IoFunctionSet
import wom.values.{WomFile, WomFloat, WomValue}

import scala.concurrent.Future
import scala.util.{Failure, Try}

case object NoIoFunctionSet extends IoFunctionSet {
  override def readFile(path: String): Future[String] = Future.failed(new NotImplementedError("readFile is not available here"))

  override def writeFile(path: String, content: String): Future[WomFile] = Future.failed(new NotImplementedError("writeFile is not available here"))

  override def stdout(params: Seq[Try[WomValue]]): Try[WomFile] = Failure(new NotImplementedError("stdout is not available here"))

  override def stderr(params: Seq[Try[WomValue]]): Try[WomFile] = Failure(new NotImplementedError("stderr is not available here"))

  override def glob(pattern: String): Future[Seq[String]] = throw new NotImplementedError("glob is not available here")

  override def size(params: Seq[Try[WomValue]]): Future[WomFloat] = Future.failed(new NotImplementedError("size is not available here"))
}
