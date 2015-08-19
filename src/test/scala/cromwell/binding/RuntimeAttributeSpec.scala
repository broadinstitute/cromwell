package cromwell.binding

import com.google.api.services.genomics.model.Disk
import cromwell.parser.{BackendType, MemorySize}
import org.scalatest.{Matchers, FlatSpec}
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

  val WorkflowWithFullGooglyConfig =
    """
      |task googly_task {
      |  command {
      |    echo "Hello JES!"
      |  }
      |  runtime {
      |    docker: "ubuntu:latest"
      |    memory: "4096 MB"
      |    cpu: "3"
      |    defaultZones: "US_Metro US_Backwater"
      |    defaultDisks: "Disk1 3 SSD, Disk2 500 OldSpinnyKind"
      |  }
      |}
      |
      |workflow googly_workflow {
      |  call googly_task
      |}
    """.stripMargin

  val WorkflowWithoutGooglyConfig =
    """
      |task googly_task {
      |  command {
      |    echo "Hello JES!"
      |  }
      |  runtime {
      |    docker: "ubuntu:latest"
      |  }
      |}
      |
      |workflow googly_workflow {
      |  call googly_task
      |}
    """.stripMargin
}

class RuntimeAttributeSpec extends FlatSpec with Matchers {
  val NamespaceWithRuntime = NamespaceWithWorkflow.load(WorkflowWithRuntime, BackendType.LOCAL)
  "WDL file with runtime" should "have runtime information" in {
    assert(NamespaceWithRuntime.workflow.calls.forall {
      _.task.runtimeAttributes.attributes.nonEmpty
    })
  }

  it should "have docker information" in {
    assert(NamespaceWithRuntime.workflow.calls forall {
      _.task.runtimeAttributes.docker.get == "ubuntu:latest"
    })
  }

  "WDL file without runtime" should "not have imported runtime information" in {
    assert(NamespaceWithWorkflow.load(WorkflowWithoutRuntime, BackendType.LOCAL).workflow.calls.forall {
      _.task.runtimeAttributes.attributes.isEmpty
    })
  }

  "WDL file with failOnStderr runtime" should "identify failOnStderr for (and only for) appropriate tasks" in {
    val NamespaceWithFailOnStderr = NamespaceWithWorkflow.load(WorkflowWithFailOnStderr, BackendType.LOCAL)
    val echoWithFailOnStderrIndex = NamespaceWithFailOnStderr.workflow.calls.indexWhere(call => call.name == "echoWithFailOnStderr")
    assert(echoWithFailOnStderrIndex >= 0)
    assert(NamespaceWithFailOnStderr.workflow.calls(echoWithFailOnStderrIndex).failOnStderr)

    val echoWithoutFailOnStderrIndex = NamespaceWithFailOnStderr.workflow.calls.indexWhere(call => call.name == "echoWithoutFailOnStderr")
    assert(echoWithoutFailOnStderrIndex >= 0)
    assert(!NamespaceWithFailOnStderr.workflow.calls(echoWithoutFailOnStderrIndex).failOnStderr)
  }

  "WDL file with Googly config" should "parse up properly" in {
    val NamespaceWithGooglyConfig = NamespaceWithWorkflow.load(WorkflowWithFullGooglyConfig, BackendType.JES)
    val calls = NamespaceWithGooglyConfig.workflow.calls
    val callIndex = calls.indexWhere(call => call.name == "googly_task")
    callIndex should be >= 0

    val googlyCall = calls(callIndex)
    val attributes = googlyCall.task.runtimeAttributes
    attributes.cpu shouldBe 3
    val firstDisk = new Disk().setName("Disk1").setSizeGb(3L).setType("SSD")
    val secondDisk = new Disk().setName("Disk2").setSizeGb(500L).setType("OldSpinnyKind")
    attributes.defaultDisks shouldEqual Seq(firstDisk, secondDisk)
    attributes.defaultZones shouldEqual Seq("US_Metro", "US_Backwater")
    attributes.memoryGB shouldBe 4
  }

  "WDL file with no Googly config" should "also parse up properly to defaults" in {
    val NamespaceWithoutGooglyConfig = NamespaceWithWorkflow.load(WorkflowWithoutGooglyConfig, BackendType.LOCAL)
    val calls = NamespaceWithoutGooglyConfig.workflow.calls
    val callIndex = calls.indexWhere(call => call.name == "googly_task")
    callIndex should be >= 0

    val googlyCall = calls(callIndex)
    val attributes = googlyCall.task.runtimeAttributes
    attributes.cpu shouldBe RuntimeAttributes.Defaults.Cpu
    attributes.defaultDisks shouldBe RuntimeAttributes.Defaults.Disk
    attributes.defaultZones shouldBe RuntimeAttributes.Defaults.Zones
    attributes.memoryGB shouldBe MemorySize.GB.fromBytes(RuntimeAttributes.Defaults.MemoryInBytes)
  }

  "WDL file with Googly config" should "issue warnings on the local backend" in {
    val workflow = NamespaceWithWorkflow.load(WorkflowWithFullGooglyConfig, BackendType.LOCAL)
    val warnings = workflow.tasks.head.runtimeAttributes.warnings
    warnings.head shouldBe "Found unsupported keys for backend 'LOCAL': cpu, defaultDisks, defaultZones, memory"
  }

  "WDL file without runtime section" should "not be accepted on JES backend as it has no docker" in {
    val ex = intercept[RuntimeException] {
      val workflow = NamespaceWithWorkflow.load(WorkflowWithoutRuntime, BackendType.JES)
    }
    ex.getMessage shouldBe "Missing required keys in runtime configuration for backend 'JES': docker"
  }

  "WDL file with runtime section but no docker" should "not be accepted on JES backend" in {
    val ex = intercept[RuntimeException] {
      val workflow = NamespaceWithWorkflow.load(WorkflowWithFailOnStderr, BackendType.JES)
    }
    ex.getMessage shouldBe "Missing required keys in runtime configuration for backend 'JES': docker"
  }
}
