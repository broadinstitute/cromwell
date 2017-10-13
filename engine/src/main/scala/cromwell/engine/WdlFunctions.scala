package cromwell.engine

import cromwell.backend.wdl.ReadLikeFunctions
import cromwell.core.path.PathBuilder
import wom.values.{WdlFile, WdlValue}

import scala.concurrent.Future
import scala.util.{Failure, Try}

class WdlFunctions(val pathBuilders: List[PathBuilder]) extends ReadLikeFunctions {
  private def fail(name: String) = Failure(new NotImplementedError(s"$name() not supported at the workflow level yet"))

  override def stdout(params: Seq[Try[WdlValue]]): Try[WdlFile] = fail("stdout")
  override def stderr(params: Seq[Try[WdlValue]]): Try[WdlFile] = fail("stderr")
  override def glob(path: String, pattern: String): Seq[String] = throw new NotImplementedError(s"glob(path, pattern) not implemented yet")

  // TODO WOM: fix this
  override def writeFile(path: String, content: String): Future[WdlFile] = Future.failed(new Exception("Can't write files"))
}
