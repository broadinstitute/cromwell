package wdl.draft2.model.examples

import wdl.draft2.model.WdlNamespaceWithWorkflow
import wom.core._

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
          throw new UnsupportedOperationException("not implemented")
      }
    }

    val ns = WdlNamespaceWithWorkflow.load(wdl, Seq(resolver _)).get

    ns.tasks foreach {task =>
      println(s"Task: ${task.name}")
    }
  }
}
