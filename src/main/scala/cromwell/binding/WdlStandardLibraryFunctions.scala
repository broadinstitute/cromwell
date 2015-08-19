package cromwell.binding

import cromwell.binding.values._

import scala.util.{Failure, Try}

trait WdlStandardLibraryFunctions extends WdlFunctions {
  private def fail(name: String) = Failure(new UnsupportedOperationException(s"$name() not implemented yet"))

  protected def stdout(params: Seq[Try[WdlValue]]): Try[WdlFile] = fail("stdout")
  protected def stderr(params: Seq[Try[WdlValue]]): Try[WdlFile] = fail("stderr")
  protected def read_lines(params: Seq[Try[WdlValue]]): Try[WdlArray] = fail("read_lines")
  protected def read_tsv(params: Seq[Try[WdlValue]]): Try[WdlArray] = fail("read_tsv")
  protected def read_map(params: Seq[Try[WdlValue]]): Try[WdlMap] = fail("read_map")
  protected def read_object(params: Seq[Try[WdlValue]]): Try[WdlObject] = fail("read_objects")
  protected def read_objects(params: Seq[Try[WdlValue]]): Try[WdlArray] = fail("read_objects")
  protected def read_json(params: Seq[Try[WdlValue]]): Try[WdlValue] = fail("read_json")
  protected def read_int(params: Seq[Try[WdlValue]]): Try[WdlInteger] = fail("read_int")
  protected def read_string(params: Seq[Try[WdlValue]]): Try[WdlString] = fail("read_string")
  protected def read_float(params: Seq[Try[WdlValue]]): Try[WdlFloat] = fail("read_float")
  protected def read_boolean(params: Seq[Try[WdlValue]]): Try[WdlBoolean] = fail("read_boolean")
  protected def write_lines(params: Seq[Try[WdlValue]]): Try[WdlFile] = fail("write_lines")
  protected def write_tsv(params: Seq[Try[WdlValue]]): Try[WdlFile] = fail("write_tsv")
  protected def write_map(params: Seq[Try[WdlValue]]): Try[WdlFile] = fail("write_map")
  protected def write_object(params: Seq[Try[WdlValue]]): Try[WdlFile] = fail("write_object")
  protected def write_objects(params: Seq[Try[WdlValue]]): Try[WdlFile] = fail("write_objects")
  protected def write_json(params: Seq[Try[WdlValue]]): Try[WdlFile] = fail("write_json")

  /**
   * Extract a single `WdlValue` from the specified `Seq`, returning `Failure` if the parameters
   * represent something other than a single `WdlValue`.
   */
  protected def extractSingleArgument(params: Seq[Try[WdlValue]]): Try[WdlValue] = {
    if (params.length != 1) Failure(new UnsupportedOperationException("Expected one argument, got " + params.length))
    else params.head
  }

  /* Returns one of the standard library functions (defined above) by name */
  def getFunction(name: String): WdlFunction = {
    val method = getClass.getMethod(name, classOf[Seq[Try[WdlValue]]])
    args => method.invoke(this, args).asInstanceOf[Try[WdlValue]]
  }
}

class NoFunctions extends WdlStandardLibraryFunctions
