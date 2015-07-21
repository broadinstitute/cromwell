package cromwell.binding

import cromwell.parser.BackendType
import org.scalatest.{FlatSpec, Matchers}

class ThreeStepImportNamespaceSpec extends FlatSpec with Matchers {
  val psTaskWdl = """
    |task ps {
    |  command {
    |    ps
    |  }
    |  output {
    |    File procs = stdout()
    |  }
    |}""".stripMargin

  val cgrepTaskWdl = """
    |task cgrep {
    |  command {
    |    grep '${pattern}' ${File in_file} | wc -l
    |  }
    |  output {
    |    Int count = read_int(stdout())
    |  }
    |}""".stripMargin

  val wcTaskWdl = """
    |task wc {
    |  command {
    |    cat ${File in_file} | wc -l
    |  }
    |  output {
    |    Int count = read_int(stdout())
    |  }
    |}""".stripMargin

  val workflowWdl = """
    |import "ps" as ns1
    |import "cgrep" as ns2
    |import "wc" as ns3
    |
    |workflow three_step {
    |  call ns1.ps
    |  call ns2.cgrep {
    |    input: in_file=ns1.ps.procs
    |  }
    |  call ns3.wc {
    |    input: in_file=ns1.ps.procs
    |  }
    |}""".stripMargin

  def resolver(importUri: String): WdlSource = {
    importUri match {
      case "ps" => psTaskWdl
      case "cgrep" => cgrepTaskWdl
      case "wc" => wcTaskWdl
      case _ => throw new RuntimeException(s"Can't resolve $importUri")
    }
  }

  val namespace = NamespaceWithWorkflow.load(workflowWdl, resolver _, BackendType.LOCAL)


  "WDL file with imports" should "Have 0 tasks (3 tasks are in separate namespace)" in {
    namespace.tasks.size shouldEqual 0
  }
  it should "Have 3 imported WdlBindings" in {
    namespace.namespaces.size shouldEqual 3
  }
  it should "Have 3 imported WdlBindings with tasks 'ps', 'cgrep', and 'wc'" in {
    namespace.namespaces flatMap {_.tasks} map {_.name} shouldEqual Seq("ps", "cgrep", "wc")
  }
  it should "Throw an exception if the import resolver fails to resolve an import" in {
    def badResolver(s: String): String = {
      throw new RuntimeException(s"Can't Resolve")
    }
    try {
      val badBinding = WdlNamespace.load(workflowWdl, badResolver _, BackendType.LOCAL)
      fail("Expecting an exception to be thrown when using badResolver")
    } catch {
      case _: RuntimeException =>
    }
  }
}

