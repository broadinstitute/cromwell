package wdl4s.examples

import wdl4s.{TaskCall, WdlNamespaceWithWorkflow}

object ex5 {
  def main(args: Array[String]) {
    val wdl = """
      |task a {
      |  command { ps }
      |  output { File procs = stdout() }
      |}
      |
      |task b {
      |  File s
      |  command { wc -l ${s} }
      |}
      |
      |workflow wf {
      | call a
      | call b {input: s=a.procs}
      |}""".stripMargin

    val ns = WdlNamespaceWithWorkflow.load(wdl, Seq.empty).get

    Seq(ns.resolve("wf.a"), ns.resolve("wf.b")) foreach {
      case Some(c: TaskCall) => println(s"Call '${c.fullyQualifiedName}' prerequisites: ${c.upstream}")
      case _ =>
    }
  }
}
