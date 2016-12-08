package wdl4s

import better.files.File
import org.scalatest.{FlatSpec, Matchers}
import wdl4s.exception.ValidationException

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
    |  String pattern
    |  String in_file
    |  command {
    |    grep '${pattern}' ${in_file} | wc -l
    |  }
    |  output {
    |    Int count = read_int(stdout())
    |  }
    |}""".stripMargin

  val wcTaskWdl = """
    |task wc {
    |  File in_file
    |  command {
    |    cat ${in_file} | wc -l
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
    |    input: in_file=ps.procs
    |  }
    |  call ns3.wc {
    |    input: in_file=ps.procs
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

  val namespace = WdlNamespaceWithWorkflow.load(workflowWdl, resolver _)


  "WDL file with imports" should "Have 0 tasks (3 tasks are in separate namespace)" in {
    namespace.tasks.size shouldEqual 0
  }
  it should "Have 3 imported WdlNamespaces" in {
    namespace.namespaces.size shouldEqual 3
  }
  it should "Have 3 imported WdlNamespaces with tasks 'ps', 'cgrep', and 'wc'" in {
    namespace.namespaces flatMap {_.tasks} map {_.name} shouldEqual Seq("ps", "cgrep", "wc")
  }
  it should "Throw an exception if the import resolver fails to resolve an import" in {
    def badResolver(s: String): String = {
      throw new RuntimeException(s"Can't Resolve")
    }
    try {
      val badBinding = WdlNamespace.loadUsingSource(workflowWdl, None, Option(Seq(badResolver)))
      fail("Expecting an exception to be thrown when using badResolver")
    } catch {
      case _: ValidationException =>
    }
  }
}

