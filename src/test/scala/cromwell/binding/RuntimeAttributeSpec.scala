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

  val WorkflowWithFailOnStderr =
    """
      |task echoWithFailOnStderr {
      |  command {
      |    echo 66555 >&2
      |  }
      |  runtime {
      |    failOnStderr: "true"
      |  }
      |}
      |task echoWithoutFailOnStderr {
      |  command {
      |    echo 66555 >&2
      |  }
      |  runtime {
      |    failOnStderr: "false"
      |  }
      |}
      |
      |workflow echo_wf {
      |  call echoWithFailOnStderr
      |  call echoWithoutFailOnStderr
      |}
    """.stripMargin

  val NamespaceWithRuntime = NamespaceWithWorkflow.load(WorkflowWithRuntime)
  val NamespaceWithoutRuntime = NamespaceWithWorkflow.load(WorkflowWithoutRuntime)
  val NamespaceWithFailOnStderr = NamespaceWithWorkflow.load(WorkflowWithFailOnStderr)
}

class RuntimeAttributeSpec extends FlatSpec with Matchers {
  "WDL file with runtime" should "have runtime information" in {
    assert(NamespaceWithRuntime.workflow.calls.forall {_.task.runtimeAttributes.nonEmpty})
  }

  it should "have docker information" in {
    assert(NamespaceWithRuntime.workflow.calls forall {_.task.runtimeAttributes.docker.get == "ubuntu:latest"})
  }

  "WDL file without runtime" should "not have imported runtime information" in {
    assert(NamespaceWithoutRuntime.workflow.calls.forall {_.task.runtimeAttributes.isEmpty})
  }

  "WDL file with failOnStderr runtime" should "identify failOnStderr for (and only for) appropriate tasks" in {
    val echoWithFailOnStderrIndex = NamespaceWithFailOnStderr.workflow.calls.indexWhere(call => call.name == "echoWithFailOnStderr")
    assert(echoWithFailOnStderrIndex >= 0)
    assert(NamespaceWithFailOnStderr.workflow.calls(echoWithFailOnStderrIndex).failOnStderr)

    val echoWithoutFailOnStderrIndex = NamespaceWithFailOnStderr.workflow.calls.indexWhere(call => call.name == "echoWithoutFailOnStderr")
    assert(echoWithoutFailOnStderrIndex >= 0)
    assert(!NamespaceWithFailOnStderr.workflow.calls(echoWithoutFailOnStderrIndex).failOnStderr)
  }
  
  
}
