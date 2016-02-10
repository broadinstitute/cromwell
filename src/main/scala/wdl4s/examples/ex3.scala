package wdl4s.examples

import wdl4s.{WdlSource, NamespaceWithWorkflow}

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

    val ns = NamespaceWithWorkflow.load(wdl, resolver)

    ns.tasks foreach {task =>
      println(s"Task: ${task.name}")
    }

    ns.namespaces foreach { n =>
      n.tasks.foreach { t =>
        println(s"Imported Task: ${t.name} (from ${n.importedAs.get})")
      }
    }
  }
}
