package cromwell.binding

import org.scalatest.{Matchers, FlatSpec}
import RuntimeAttributeSpec._

object RuntimeAttributeSpec {
  val WorkflowWithRuntime = """
                              |task ps {
                              |  command {
                              |    ps
                              |  }
                              |  output {
                              |    File procs = stdout()
                              |  }
                              |  runtime {
                              |    docker: "ubuntu:latest"
                              |  }
                              |}
                              |
                              |task cgrep {
                              |  command {
                              |    grep '${pattern}' ${File in_file} | wc -l
                              |  }
                              |  output {
                              |    Int count = read_int(stdout())
                              |  }
                              |  runtime {
                              |    docker: "ubuntu:latest"
                              |  }
                              |}
                              |
                              |task wc {
                              |  command {
                              |    cat ${File in_file} | wc -l
                              |  }
                              |  output {
                              |    Int count = read_int(stdout())
                              |  }
                              |  runtime {
                              |     docker: "ubuntu:latest"
                              |  }
                              |}
                              |
                              |workflow three_step {
                              |  call ps
                              |  call cgrep {
                              |    input: in_file=ps.procs
                              |  }
                              |  call wc {
                              |    input: in_file=ps.procs
                              |  }
                              |}
                              |
                            """.stripMargin


  val WorkflowWithoutRuntime = """
                                 |task hello {
                                 |  command {
                                 |    echo "Hello ${addressee}!"
                                 |  }
                                 |  output {
                                 |    String salutation = read_string(stdout())
                                 |  }
                                 |}
                                 |
                                 |workflow hello {
                                 |  call hello
                                 |}
                               """.stripMargin

  val NamespaceWithRuntime = WdlNamespace.load(WorkflowWithRuntime)
  val NamespaceWithoutRuntime = WdlNamespace.load(WorkflowWithoutRuntime)
}

class RuntimeAttributeSpec extends FlatSpec with Matchers {
  "WDL file with runtime" should "have runtime information" in {
    assert(NamespaceWithRuntime.workflows.head.calls.forall {_.task.runtimeAttributes.nonEmpty})
  }

  it should "have docker information" in {
    assert(NamespaceWithRuntime.workflows.head.calls forall {_.task.runtimeAttributes.docker.get == "ubuntu:latest"})
  }

  "WDL file without runtime" should "not have imported runtime information" in {
    assert(NamespaceWithoutRuntime.workflows.head.calls.forall {_.task.runtimeAttributes.isEmpty})
  }
}
