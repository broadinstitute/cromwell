package wdl4s.examples

import wdl4s.NamespaceWithWorkflow
import wdl4s.expression.NoFunctions
import wdl4s.values.{WdlString, WdlValue}

object ex6 {
  def main(args: Array[String]) {
    val wdl = """
      |workflow wf {
      |  String a = "foo" + "bar"
      |  String b = "hello " + variable
      |  String c = "hello " + other_variable
      |}""".stripMargin

    val ns = NamespaceWithWorkflow.load(wdl)
    def lookup(name: String): WdlValue = {
      name match {
        case "variable" => WdlString("world")
        case _ => throw new NoSuchElementException
      }
    }
    ns.workflow.declarations foreach { decl =>
      val value = decl.expression.get.evaluate(lookup, NoFunctions)
      println(s"Declaration '${decl.toWdlString}' evaluates to: ${value}")
    }
  }
}
