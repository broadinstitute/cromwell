package cromwell.backend.wdl

import wdl4s.expression.WdlStandardLibraryFunctions
import wdl4s.values._

import scala.language.{implicitConversions, postfixOps}
import scala.util.{Failure, Try}

@deprecated("This trait will not be used in PBE world")
trait WdlStandardLibraryImpl extends WdlStandardLibraryFunctions with ReadLikeFunctions with WriteFunctions with PureFunctions {
  private def fail(name: String) = Failure(new NotImplementedError(s"$name() not implemented yet"))

  override def stdout(params: Seq[Try[WdlValue]]): Try[WdlFile] = fail("stdout")
  override def stderr(params: Seq[Try[WdlValue]]): Try[WdlFile] = fail("stderr")
}
