package wdl4s.wdl.examples

import wdl4s.wdl.{WdlTaskCall, WdlNamespaceWithWorkflow}

object ex5 {
  def main(args: Array[String]): Unit = {
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
      case Some(c: WdlTaskCall) => println(s"Call '${c.fullyQualifiedName}' prerequisites: ${c.upstream}")
      case _ =>
    }
  }
}
