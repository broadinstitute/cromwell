package wdl4s

import org.scalatest.{EitherValues, FlatSpec, Matchers}

object RuntimeAttributeSpec {
  val WorkflowWithRuntime =
    """
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
      |  String pattern
      |  File in_file
      |  command {
      |    grep '${pattern}' ${in_file} | wc -l
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
      |  File in_file
      |  command {
      |    cat ${in_file} | wc -l
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

  val WorkflowWithoutRuntime =
    """
      |task hello {
      |  String addressee
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

  class RuntimeAttributeSpec extends FlatSpec with Matchers with EitherValues {
    val NamespaceWithRuntime = WdlNamespaceWithWorkflow.load(WorkflowWithRuntime)
    val NamespaceWithoutRuntime = WdlNamespaceWithWorkflow.load(WorkflowWithoutRuntime)

    "WDL file with runtime attributes" should "have attribute maps" in {
      NamespaceWithRuntime.tasks.forall(_.runtimeAttributes.attrs.nonEmpty) should be(true)
    }

    "WDL file without runtime attributes" should "not have attribute maps" in {
      NamespaceWithoutRuntime.tasks.forall(_.runtimeAttributes.attrs.isEmpty) should be(true)
    }
  }
}
