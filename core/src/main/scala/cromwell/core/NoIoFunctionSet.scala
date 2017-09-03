package cromwell.core

import wdl4s.wdl.values.{WdlFile, WdlFloat, WdlValue}
import wdl4s.wom.expression.IoFunctionSet

import scala.concurrent.Future
import scala.util.Try

case object NoIoFunctionSet extends IoFunctionSet {
  override def readFile(path: String): Future[String] = ???

  override def writeFile(path: String, content: String): Future[WdlFile] = ???

  override def stdout(params: Seq[Try[WdlValue]]): Try[WdlFile] = ???

  override def stderr(params: Seq[Try[WdlValue]]): Try[WdlFile] = ???

  override def glob(path: String, pattern: String): Seq[String] = ???

  override def size(params: Seq[Try[WdlValue]]): Try[WdlFloat] = ???
}