package wdl.draft2.model.examples

import wdl.draft2.Draft2ResolvedImportBundle
import wdl.draft2.model.WdlNamespaceWithWorkflow
import wom.ResolvedImportRecord

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

    def resolver(importString: String): Draft2ResolvedImportBundle = {
      importString match {
        case "some_string" => Draft2ResolvedImportBundle("task imported { command {ps} }", ResolvedImportRecord("some_string"))
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
