package cromwell.engine

import java.nio.file.FileSystem

import cromwell.backend.wdl.{PureFunctions, ReadLikeFunctions}
import wdl4s.expression.WdlStandardLibraryFunctions
import wdl4s.values.{WdlFile, WdlValue}

import scala.util.{Failure, Try}

class WdlFunctions(val fileSystems: List[FileSystem]) extends WdlStandardLibraryFunctions with ReadLikeFunctions with PureFunctions {
  /**
    * Ordered list of filesystems to be used to execute WDL functions needing IO.
    */
  private def fail(name: String) = Failure(new NotImplementedError(s"$name() not supported at the workflow level yet"))

  override def write_json(params: Seq[Try[WdlValue]]): Try[WdlFile] = fail("write_json")
  override def stdout(params: Seq[Try[WdlValue]]): Try[WdlFile] = fail("stdout")
  override def stderr(params: Seq[Try[WdlValue]]): Try[WdlFile] = fail("stderr")
  override def glob(path: String, pattern: String): Seq[String] = throw new NotImplementedError(s"glob(path, pattern) not implemented yet")
  override def writeTempFile(path: String, prefix: String, suffix: String, content: String): String = throw new NotImplementedError(s"Cromwell doesn't support write_* functions at the workflow level")
  override def write_tsv(params: Seq[Try[WdlValue]]): Try[WdlFile] = fail("write_tsv")
}
