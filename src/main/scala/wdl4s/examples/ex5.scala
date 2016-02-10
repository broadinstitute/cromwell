package wdl4s.examples

import wdl4s.{Call, NamespaceWithWorkflow}

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

    val ns = NamespaceWithWorkflow.load(wdl)

    Seq(ns.resolve("wf.a"), ns.resolve("wf.b")) foreach { scope =>
      scope match {
        case Some(c: Call) => println(s"Call '${c.fullyQualifiedName}' prerequisites: ${c.prerequisiteScopes}")
      }
    }
  }
}
