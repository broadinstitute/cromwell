package wdl4s.wdl

import org.scalatest.{FlatSpec, Matchers}
import wdl4s.wdl.exception.ValidationException

import scala.util.Failure

class ThreeStepImportSpec extends FlatSpec with Matchers {
  val psTaskWdl = """
    |task ps {
    |  command {
    |    ps
    |  }
    |  output {
    |    File procs = stdout()
    |  }
    |}""".stripMargin

  val cgrepTaskWdl = s"""
    |task cgrep {
    |  String pattern
    |  File in_file
    |  command {
    |    grep '$${pattern}' $${in_file} | wc -l
    |  }
    |  output {
    |    Int count = read_int(stdout())
    |  }
    |}""".stripMargin

  val wcTaskWdl = s"""
    |task wc {
    |  File in_file
    |  command {
    |    cat $${in_file} | wc -l
    |  }
    |  output {
    |    Int count = read_int(stdout())
    |  }
    |}""".stripMargin

  val workflowWdl = """
    |import "ps"
    |import "cgrep"
    |import "wc"
    |
    |workflow three_step {
    |  call ps.ps
    |  call cgrep.cgrep {
    |    input: in_file=ps.procs
    |  }
    |  call wc.wc {
    |    input: in_file=ps.procs
    |  }
    |}""".stripMargin

  def resolver(importUri: String): String = {
    importUri match {
      case "ps" => psTaskWdl
      case "cgrep" => cgrepTaskWdl
      case "wc" => wcTaskWdl
      case _ => throw new RuntimeException(s"Can't resolve $importUri")
    }
  }

  val namespace = WdlNamespaceWithWorkflow.load(workflowWdl, Seq(resolver _)).get

  "WDL file with imports" should "Have 3 tasks" in {
    namespace.tasks.size shouldEqual 0
  }

  it should "Have 0 imported WdlNamespaces" in {
    namespace.namespaces.size shouldEqual 3
  }
  it should "Have imported namespaces with tasks named 'ps', 'cgrep' and 'wc'" in {
    namespace.namespaces flatMap {_.tasks} map {_.name} shouldEqual Seq("ps", "cgrep", "wc")
  }
  it should "Throw an exception if the import resolver fails to resolve an import" in {
    def badResolver(s: String): String = {
      throw new RuntimeException(s"Can't Resolve")
    }

    WdlNamespace.loadUsingSource(workflowWdl, None, Option(Seq(badResolver))) match {
      case Failure(_: ValidationException) =>
      case x => fail(s"Expecting an exception to be thrown when using badResolver but got $x")
    }
  }
}
