package wdl4s.wdl.examples
import wdl4s.wdl.WdlNamespaceWithWorkflow

object ex1 {
  def main(args: Array[String]): Unit = {
    val wdl = """
    |task a {
    |  command { ps }
    |}
    |workflow wf {
    | call a
    |}""".stripMargin

    val ns = WdlNamespaceWithWorkflow.load(wdl, Seq.empty).get

    println(s"Workflow: ${ns.workflow.unqualifiedName}")
    ns.workflow.calls foreach {call =>
      println(s"Call: ${call.unqualifiedName}")
    }

    ns.tasks foreach {task =>
      println(s"Task: ${task.name}")
      println(s"Command: ${task.commandTemplate}")
    }
  }
}
