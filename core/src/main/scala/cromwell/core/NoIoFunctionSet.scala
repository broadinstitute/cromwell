package cromwell.core

import wom.expression.IoFunctionSet
import wom.values.{WdlFile, WdlFloat, WdlValue}

import scala.concurrent.Future
import scala.util.{Failure, Try}

case object NoIoFunctionSet extends IoFunctionSet {
  override def readFile(path: String): Future[String] = Future.failed(new NotImplementedError("readFile is not available here"))

  override def writeFile(path: String, content: String): Future[WdlFile] = Future.failed(new NotImplementedError("writeFile is not available here"))

  override def stdout(params: Seq[Try[WdlValue]]): Try[WdlFile] = Failure(new NotImplementedError("stdout is not available here"))

  override def stderr(params: Seq[Try[WdlValue]]): Try[WdlFile] = Failure(new NotImplementedError("stderr is not available here"))

  override def glob(path: String, pattern: String): Seq[String] = throw new NotImplementedError("glob is not available here")

  override def size(params: Seq[Try[WdlValue]]): Try[WdlFloat] = Failure(new NotImplementedError("size is not available here"))
}
