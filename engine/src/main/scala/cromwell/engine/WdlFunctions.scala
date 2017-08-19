package cromwell.engine

import cromwell.backend.wdl.ReadLikeFunctions
import wdl4s.wdl.expression.PureStandardLibraryFunctionsLike
import cromwell.core.path.PathBuilder
import wdl4s.wdl.values.{WdlFile, WdlValue}

import scala.util.{Failure, Try}

class WdlFunctions(val pathBuilders: List[PathBuilder]) extends PureStandardLibraryFunctionsLike with ReadLikeFunctions {
  private def fail(name: String) = Failure(new NotImplementedError(s"$name() not supported at the workflow level yet"))

  override def write_json(params: Seq[Try[WdlValue]]): Try[WdlFile] = fail("write_json")
  override def stdout(params: Seq[Try[WdlValue]]): Try[WdlFile] = fail("stdout")
  override def stderr(params: Seq[Try[WdlValue]]): Try[WdlFile] = fail("stderr")
  override def glob(path: String, pattern: String): Seq[String] = throw new NotImplementedError(s"glob(path, pattern) not implemented yet")
  override def write_tsv(params: Seq[Try[WdlValue]]): Try[WdlFile] = fail("write_tsv")
}
