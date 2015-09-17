package cromwell.binding.expression

import cromwell.binding.WdlExpressionException
import cromwell.binding.types._
import cromwell.binding.values._

import scala.util.{Failure, Success, Try}

trait WdlStandardLibraryFunctions extends WdlFunctions[WdlValue] {
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
}

class NoFunctions extends WdlStandardLibraryFunctions

class WdlStandardLibraryFunctionsType extends WdlFunctions[WdlType] {
  protected def stdout(params: Seq[Try[WdlType]]): Try[WdlType] = Success(WdlFileType)
  protected def stderr(params: Seq[Try[WdlType]]): Try[WdlType] = Success(WdlFileType)
  protected def read_lines(params: Seq[Try[WdlType]]): Try[WdlType] = Success(WdlArrayType(WdlStringType))
  protected def read_tsv(params: Seq[Try[WdlType]]): Try[WdlType] = Success(WdlArrayType(WdlArrayType(WdlStringType)))
  protected def read_map(params: Seq[Try[WdlType]]): Try[WdlType] = Success(WdlMapType(WdlStringType, WdlStringType))
  protected def read_object(params: Seq[Try[WdlType]]): Try[WdlType] = Success(WdlObjectType)
  protected def read_objects(params: Seq[Try[WdlType]]): Try[WdlType] = Success(WdlArrayType(WdlObjectType))
  protected def read_json(params: Seq[Try[WdlType]]): Try[WdlType] = Failure(new WdlExpressionException("Return type of read_json() can't be known statically"))
  protected def read_int(params: Seq[Try[WdlType]]): Try[WdlType] = Success(WdlIntegerType)
  protected def read_string(params: Seq[Try[WdlType]]): Try[WdlType] = Success(WdlStringType)
  protected def read_float(params: Seq[Try[WdlType]]): Try[WdlType] = Success(WdlFloatType)
  protected def read_boolean(params: Seq[Try[WdlType]]): Try[WdlType] = Success(WdlBooleanType)
  protected def write_lines(params: Seq[Try[WdlType]]): Try[WdlType] = Success(WdlFileType)
  protected def write_tsv(params: Seq[Try[WdlType]]): Try[WdlType] = Success(WdlFileType)
  protected def write_map(params: Seq[Try[WdlType]]): Try[WdlType] = Success(WdlFileType)
  protected def write_object(params: Seq[Try[WdlType]]): Try[WdlType] = Success(WdlFileType)
  protected def write_objects(params: Seq[Try[WdlType]]): Try[WdlType] = Success(WdlFileType)
  protected def write_json(params: Seq[Try[WdlType]]): Try[WdlType] = Success(WdlFileType)
}
