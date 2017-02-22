package wdl4s.examples

import wdl4s.{WdlSource, WdlNamespaceWithWorkflow}

object ex3 {
  def main(args: Array[String]) {
    val wdl = """
      |import "some_string" as my_namespace
      |task a {
      |  command { ps }
      |}
      |workflow wf {
      | call a
      |}""".stripMargin

    def resolver(importString: String): WdlSource = {
      importString match {
        case "some_string" => "task imported { command {ps} }"
        case _ => throw new NotImplementedError()
      }
    }

    val ns = WdlNamespaceWithWorkflow.load(wdl, Seq(resolver _)).get

    ns.tasks foreach {task =>
      println(s"Task: ${task.name}")
    }

    ns.namespaces foreach { n =>
      n.tasks.foreach { t =>
        println(s"Imported Task: ${t.name} (from ${n.importedAs})")
      }
    }
  }
}
