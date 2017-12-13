package wdl.examples

import common.validation.Validation._
import wdl.WdlNamespaceWithWorkflow
import wdl.expression.WdlFunctions
import wom.types.{WomArrayType, WomIntegerType}
import wom.values.{WomArray, WomInteger, WomSingleFile, WomString, WomValue}

import scala.util.{Success, Try}

object ex7 {
  def main(args: Array[String]): Unit = {
    val wdl = s"""
      |task a {
      |  String prefix
      |  Array[Int] ints
      |  command {
      |    python script.py $${write_lines(ints)} > $${prefix + ".out"}
      |  }
      |}
      |workflow wf {
      |  call a
      |}""".stripMargin

    val ns = WdlNamespaceWithWorkflow.load(wdl, Seq.empty).get
    val inputs = Map(
      "prefix" -> WomString("some_prefix"),
      "ints" -> WomArray(WomArrayType(WomIntegerType), Seq(1,2,3,4,5).map(WomInteger.apply))
    )

    class CustomFunctions extends WdlFunctions[WomValue] {
      def write_lines(params: Seq[Try[WomValue]]): Try[WomValue] = {
        // Validate `params`, write the result to a file, return file path
        Success(WomSingleFile("/tmp/array.txt"))
      }
    }

    ns.taskCalls.find( _.unqualifiedName == "a") foreach { call =>
      val wdlFunctions: CustomFunctions = new CustomFunctions
      val evaluatedInputs = call.evaluateTaskInputs(inputs, wdlFunctions).get
      println(call.task.instantiateCommand(evaluatedInputs, wdlFunctions).toTry.get)
    }
  }
}
