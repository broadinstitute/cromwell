package wdl.examples

import wdl.WdlNamespaceWithWorkflow
import wdl.expression.NoFunctions
import wom.values.{WomString, WomValue}

object ex6 {
  def main(args: Array[String]): Unit = {
    val wdl = """
      |workflow wf {
      |  String a = "foo" + "bar"
      |  String b = "hello " + variable
      |  String c = "hello " + other_variable
      |}""".stripMargin

    val ns = WdlNamespaceWithWorkflow.load(wdl, Seq.empty).get
    def lookup(name: String): WomValue = {
      name match {
        case "variable" => WomString("world")
        case _ => throw new NoSuchElementException
      }
    }
    ns.workflow.declarations foreach { decl =>
      val value = decl.expression.get.evaluate(lookup, NoFunctions)
      println(s"Declaration '${decl.toWdlString}' evaluates to: $value")
    }
  }
}
