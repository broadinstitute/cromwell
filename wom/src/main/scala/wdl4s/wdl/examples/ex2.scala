package wdl4s.wdl.examples

import wdl4s.wdl.{WorkflowSource, WdlNamespaceWithWorkflow}

object ex2 {
  def main(args: Array[String]): Unit = {
    val wdl = """
      |import "some_string"
      |task a {
      |  command { ps }
      |}
      |workflow wf {
      | call a
      |}""".stripMargin

    def resolver(importString: String): WorkflowSource = {
      importString match {
        case "some_string" => "task imported { command {ps} }"
        case s if s.startsWith("http://") =>
          // issue HTTP request
          throw new NotImplementedError("not implemented")
      }
    }

    val ns = WdlNamespaceWithWorkflow.load(wdl, Seq(resolver _)).get

    ns.tasks foreach {task =>
      println(s"Task: ${task.name}")
    }
  }
}
