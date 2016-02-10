package wdl4s.examples

import wdl4s.{WdlSource, NamespaceWithWorkflow}

object ex2 {
  def main(args: Array[String]) {
    val wdl = """
      |import "some_string"
      |task a {
      |  command { ps }
      |}
      |workflow wf {
      | call a
      |}""".stripMargin

    def resolver(importString: String): WdlSource = {
      importString match {
        case "some_string" => "task imported { command {ps} }"
        case s if s.startsWith("http://") =>
          // issue HTTP request
          throw new NotImplementedError("not implemented")
      }
    }

    val ns = NamespaceWithWorkflow.load(wdl, resolver)

    ns.tasks foreach {task =>
      println(s"Task: ${task.name}")
    }
  }
}
