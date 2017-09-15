package cromwell.core

import wdl4s.wdl.values.{WdlFile, WdlFloat, WdlValue}
import wdl4s.wom.expression.IoFunctionSet

import scala.concurrent.Future
import scala.util.Try

case object NoIoFunctionSet extends IoFunctionSet {
  override def readFile(path: String): Future[String] = throw new NotImplementedError("readFile is not available here")

  override def writeFile(path: String, content: String): Future[WdlFile] = throw new NotImplementedError("writeFile is not available here")

  override def stdout(params: Seq[Try[WdlValue]]): Try[WdlFile] = throw new NotImplementedError("stdout is not available here")

  override def stderr(params: Seq[Try[WdlValue]]): Try[WdlFile] = throw new NotImplementedError("stderr is not available here")

  override def glob(path: String, pattern: String): Seq[String] = throw new NotImplementedError("glob is not available here")

  override def size(params: Seq[Try[WdlValue]]): Try[WdlFloat] = throw new NotImplementedError("size is not available here")
}
