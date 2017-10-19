package wdl.examples

import wdl.WdlNamespaceWithWorkflow
import wdl.expression.WdlFunctions
import wom.types.{WdlArrayType, WdlIntegerType}
import wom.values.{WdlArray, WdlFile, WdlInteger, WdlString, WdlValue}

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
      "prefix" -> WdlString("some_prefix"),
      "ints" -> WdlArray(WdlArrayType(WdlIntegerType), Seq(1,2,3,4,5).map(WdlInteger(_)))
    )

    class CustomFunctions extends WdlFunctions[WdlValue] {
      def write_lines(params: Seq[Try[WdlValue]]): Try[WdlValue] = {
        // Validate `params`, write the result to a file, return file path
        Success(WdlFile("/tmp/array.txt"))
      }
    }

    ns.taskCalls.find( _.unqualifiedName == "a") foreach { call =>
      val wdlFunctions: CustomFunctions = new CustomFunctions
      val evaluatedInputs = call.evaluateTaskInputs(inputs, wdlFunctions).get
      println(call.task.instantiateCommand(evaluatedInputs, wdlFunctions).get)
    }
  }
}
