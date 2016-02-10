package cromwell.engine

import wdl4s.expression.NoFunctions
import wdl4s.values.WdlValue
import org.scalatest.{FlatSpec, Matchers}

import scala.util.{Failure, Success, Try}

class EngineFunctionsSpec extends FlatSpec with Matchers {
  def expectFailure(value: Try[WdlValue]) = value match {
    case Success(s) => fail(s"$s: Expected this function invocation to fail")
    case Failure(ex) => // expected
  }
  "EngineFunctions" should "all initially be undefined" in {
    val stdFunctions = Seq(
      "stdout", "stderr", "read_lines", "read_tsv", "read_map", "read_object", "read_objects",
      "read_json", "read_int", "read_string", "read_float", "read_boolean", "write_lines",
      "write_tsv", "write_map", "write_object", "write_objects", "write_json", "size"
    )
    stdFunctions.foreach {func =>
      expectFailure(NoFunctions.getFunction(func)(Seq.empty[Try[WdlValue]]))
    }
  }
}
