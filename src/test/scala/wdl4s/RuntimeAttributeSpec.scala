package wdl4s

import org.scalatest.{EitherValues, Matchers, FlatSpec}
import RuntimeAttributeSpec._

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

  val WorkflowWithMessedUpMemory =
    """
    |task messed_up_memory {
    |  command {
    |      echo "YO"
    |  }
    |  runtime {
    |    memory: "HI TY"
    |  }
    |}
    |
    |workflow great_googly_moogly {
    |  call messed_up_memory
    |}
  """.stripMargin
  val WorkflowWithMessedUpMemoryUnit =
    """
      |task messed_up_memory {
      |  command {
      |      echo "YO"
      |  }
      |  runtime {
      |    memory: "5 TY"
      |  }
      |}
      |
      |workflow great_googly_moogly {
      |  call messed_up_memory
      |}
    """.stripMargin
  }

class RuntimeAttributeSpec extends FlatSpec with Matchers with EitherValues {
  val NamespaceWithRuntime = NamespaceWithWorkflow.load(WorkflowWithRuntime)
  val NamespaceWithoutRuntime = NamespaceWithWorkflow.load(WorkflowWithoutRuntime)

  "WDL file with runtime attributes" should "have attribute maps" in {
    NamespaceWithRuntime.tasks.forall(_.runtimeAttributes.attrs.nonEmpty) should be(true)
  }

  "WDL file without runtime attributes" should "not have attribute maps" in {
    NamespaceWithoutRuntime.tasks.forall(_.runtimeAttributes.attrs.isEmpty) should be(true)
  }

  "WDL file with a seriously screwed up memory runtime" should "not parse" in {
    val ex = intercept[IllegalArgumentException] {
      val namespaceWithBorkedMemory = NamespaceWithWorkflow.load(WorkflowWithMessedUpMemory)
    }

    ex.getMessage should include ("should be of the form X Unit")
  }

  "WDL file with an invalid memory unit" should "say so" in {
    val ex = intercept[IllegalArgumentException] {
      val namespaceWithBorkedMemory = NamespaceWithWorkflow.load(WorkflowWithMessedUpMemoryUnit)
    }

    ex.getMessage should include ("is an invalid memory unit")
  }
}
