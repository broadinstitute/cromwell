package wdl4s.wdl.examples

import wdl4s.wdl.WdlNamespaceWithWorkflow

object ex4 {
  def main(args: Array[String]): Unit = {
    val wdl = """
      |task a {
      |  command { ps }
      |}
      |workflow wf {
      | call a
      | call a as b
      |}""".stripMargin

    val ns = WdlNamespaceWithWorkflow.load(wdl, Seq.empty).get

    println(ns.resolve("wf.a")) // resolves to Call object for `call a`
    println(ns.resolve("wf.b")) // resolves to Call object for `call a as b`
    println(ns.findTask("a")) // resolves to Task object for `task a`
  }
}
