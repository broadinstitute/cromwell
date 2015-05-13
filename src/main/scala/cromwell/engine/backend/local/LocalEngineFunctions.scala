package cromwell.engine.backend.local

import cromwell.binding.values.{WdlFile, WdlInteger, WdlString, WdlValue}
import cromwell.engine.EngineFunctions

import scala.util.{Success, Try, Failure}

class LocalEngineFunctions(executionContext: TaskExecutionContext) extends EngineFunctions {

  def assertSingleArgument(params: Seq[Try[WdlValue]]): Try[WdlValue] = {
    if (params.length != 1) Failure(new UnsupportedOperationException("Expected one argument, got " + params.length))
    else params.head
  }

  override protected def read_int(params: Seq[Try[WdlValue]]): Try[WdlInteger] = {
    assertSingleArgument(params).map {
      case f: WdlFile =>
        WdlInteger(io.Source.fromFile(f.value).mkString.toInt)
      case e =>
        throw new UnsupportedOperationException("Unsupported value type " + e)
    }
  }

  override protected def read_string(params: Seq[Try[WdlValue]]): Try[WdlString] = {
    assertSingleArgument(params).map {
      case f: WdlFile =>
        WdlString(io.Source.fromFile(f.value).mkString)
      case s: WdlString if s.value == "stdout" =>
        read_string(Seq(Success(WdlFile(executionContext.stdout)))).get
      case e =>
        throw new UnsupportedOperationException("Unsupported argument " + e)
    }
  }
}
