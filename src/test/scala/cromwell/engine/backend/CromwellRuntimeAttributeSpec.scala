package cromwell.engine.backend

import com.google.api.services.genomics.model.Disk
import cromwell.engine.backend.runtimeattributes.{AttributeMap, CromwellRuntimeAttributes, ContinueOnReturnCodeFlag, ContinueOnReturnCodeSet}
import wdl4s.NamespaceWithWorkflow
import cromwell.engine.backend.CromwellRuntimeAttributeSpec._
import org.scalatest.{EitherValues, FlatSpec, Matchers}

object CromwellRuntimeAttributeSpec {
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

  val WorkflowWithContinueOnReturnCode =
    """
      |task echoWithSingleContinueOnReturnCode {
      |  command {
      |    cat non_existent_file
      |  }
      |  runtime {
      |    continueOnReturnCode: 123
      |  }
      |}
      |task echoWithExpressionContinueOnReturnCode {
      |  command {
      |    cat non_existent_file
      |  }
      |  runtime {
      |    continueOnReturnCode: 123 + 321
      |  }
      |}
      |task echoWithListContinueOnReturnCode {
      |  command {
      |    cat non_existent_file
      |  }
      |  runtime {
      |    continueOnReturnCode: [0, 1, 2, 3]
      |  }
      |}
      |task echoWithTrueContinueOnReturnCode {
      |  command {
      |    cat non_existent_file
      |  }
      |  runtime {
      |    continueOnReturnCode: true
      |  }
      |}
      |task echoWithFalseContinueOnReturnCode {
      |  command {
      |    cat non_existent_file
      |  }
      |  runtime {
      |    continueOnReturnCode: false
      |  }
      |}
      |task echoWithTrueStringContinueOnReturnCode {
      |  command {
      |    cat non_existent_file
      |  }
      |  runtime {
      |    continueOnReturnCode: "true"
      |  }
      |}
      |task echoWithFalseStringContinueOnReturnCode {
      |  command {
      |    cat non_existent_file
      |  }
      |  runtime {
      |    continueOnReturnCode: "false"
      |  }
      |}
      |
      |workflow echo_wf {
      |  call echoWithSingleContinueOnReturnCode
      |  call echoWithExpressionContinueOnReturnCode
      |  call echoWithListContinueOnReturnCode
      |  call echoWithTrueContinueOnReturnCode
      |  call echoWithFalseContinueOnReturnCode
      |  call echoWithTrueStringContinueOnReturnCode
      |  call echoWithFalseStringContinueOnReturnCode
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
      |    memory: "4G"
      |    cpu: "3"
      |    defaultZones: "US_Metro US_Backwater"
      |    defaultDisks: "Disk1 3 SSD, Disk2 500 HDD"
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

  val WorkflowWithLocalDiskGooglyConfig =
    """
      |task googly_task {
      |  command {
      |    echo "Hello JES!"
      |  }
      |  runtime {
      |    docker: "ubuntu:latest"
      |    defaultDisks: "local-disk 123 HDD"
      |  }
      |}
      |
      |workflow googly_workflow {
      |  call googly_task
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

  val WorkflowWithMessedUpLocalDisk =
    """
      |task messed_up_disk {
      |  command {
      |      echo "YO"
      |  }
      |  runtime {
      |    docker: "ubuntu:latest"
      |    defaultDisks: "Disk1 123 LOCAL"
      |  }
      |}
      |
      |workflow great_googly_moogly {
      |  call messed_up_disk
      |}
    """.stripMargin

  val WorkflowWithMessedUpDiskSize =
    """
      |task messed_up_disk {
      |  command {
      |      echo "YO"
      |  }
      |  runtime {
      |    docker: "ubuntu:latest"
      |    defaultDisks: "Disk1 123.0 SSD"
      |  }
      |}
      |
      |workflow great_googly_moogly {
      |  call messed_up_disk
      |}
    """.stripMargin

  val WorkflowWithMessedUpDiskType =
    """
      |task messed_up_disk {
      |  command {
      |      echo "YO"
      |  }
      |  runtime {
      |    docker: "ubuntu:latest"
      |    defaultDisks: "Disk1 123 SDD"
      |  }
      |}
      |
      |workflow great_googly_moogly {
      |  call messed_up_disk
      |}
    """.stripMargin

  implicit class EnhancedNamespace(val namespace: NamespaceWithWorkflow) extends AnyVal {
    def toCromwellRuntimeAttribute(index: Int): CromwellRuntimeAttributes = toCromwellRuntimeAttribute(index, BackendType.LOCAL)

    def toCromwellRuntimeAttribute(index: Int, backendType: BackendType): CromwellRuntimeAttributes = {
      CromwellRuntimeAttributes(namespace.workflow.calls(index).task.runtimeAttributes, backendType)
    }

    def toCromwellRuntimeAttributes(backendType: BackendType) = {
      namespace.workflow.calls map { x => CromwellRuntimeAttributes(x.task.runtimeAttributes, backendType) }
    }
  }
}

class CromwellRuntimeAttributeSpec extends FlatSpec with Matchers with EitherValues {
  val NamespaceWithRuntime = NamespaceWithWorkflow.load(WorkflowWithRuntime)

  it should "have docker information" in {
    assert(NamespaceWithRuntime.workflow.calls forall { x =>
      CromwellRuntimeAttributes(x.task.runtimeAttributes, BackendType.LOCAL).docker.get == "ubuntu:latest"
    })
  }

  "WDL file with failOnStderr runtime" should "identify failOnStderr for (and only for) appropriate tasks" in {
    val namespaceWithFailOnStderr = NamespaceWithWorkflow.load(WorkflowWithFailOnStderr)
    val echoWithFailOnStderrIndex = namespaceWithFailOnStderr.workflow.calls.indexWhere(call => call.unqualifiedName == "echoWithFailOnStderr")
    assert(echoWithFailOnStderrIndex >= 0)
    assert(namespaceWithFailOnStderr.toCromwellRuntimeAttribute(echoWithFailOnStderrIndex).failOnStderr)

    val echoWithoutFailOnStderrIndex = namespaceWithFailOnStderr.workflow.calls.indexWhere(call => call.unqualifiedName == "echoWithoutFailOnStderr")
    assert(echoWithoutFailOnStderrIndex >= 0)
    assert(!namespaceWithFailOnStderr.toCromwellRuntimeAttribute(echoWithoutFailOnStderrIndex).failOnStderr)
  }

  "WDL file with continueOnReturnCode runtime" should "identify continueOnReturnCode for (and only for) appropriate tasks" in {
    val namespaceWithContinueOnReturnCode = NamespaceWithWorkflow.load(WorkflowWithContinueOnReturnCode)

    val echoWithSingleContinueOnReturnCodeIndex =
      namespaceWithContinueOnReturnCode.workflow.calls indexWhere { call =>
        call.unqualifiedName == "echoWithSingleContinueOnReturnCode"
      }
    echoWithSingleContinueOnReturnCodeIndex should be >= 0
    namespaceWithContinueOnReturnCode.toCromwellRuntimeAttribute(echoWithSingleContinueOnReturnCodeIndex).continueOnReturnCode should be (ContinueOnReturnCodeSet(Set(123)))

    val echoWithExpressionContinueOnReturnCodeIndex =
      namespaceWithContinueOnReturnCode.workflow.calls indexWhere { call =>
        call.unqualifiedName == "echoWithExpressionContinueOnReturnCode"
      }
    echoWithExpressionContinueOnReturnCodeIndex should be >= 0
    namespaceWithContinueOnReturnCode.toCromwellRuntimeAttribute(echoWithExpressionContinueOnReturnCodeIndex).continueOnReturnCode should be(ContinueOnReturnCodeSet(Set(444)))

    val echoWithListContinueOnReturnCodeIndex =
      namespaceWithContinueOnReturnCode.workflow.calls indexWhere { call =>
        call.unqualifiedName == "echoWithListContinueOnReturnCode"
      }
    echoWithListContinueOnReturnCodeIndex should be >= 0
    namespaceWithContinueOnReturnCode.toCromwellRuntimeAttribute(echoWithListContinueOnReturnCodeIndex).continueOnReturnCode should be(ContinueOnReturnCodeSet(Set(0, 1, 2, 3)))

    val echoWithTrueContinueOnReturnCodeIndex =
      namespaceWithContinueOnReturnCode.workflow.calls indexWhere { call =>
        call.unqualifiedName == "echoWithTrueContinueOnReturnCode"
      }
    echoWithTrueContinueOnReturnCodeIndex should be >= 0
    namespaceWithContinueOnReturnCode.toCromwellRuntimeAttribute(echoWithTrueContinueOnReturnCodeIndex).continueOnReturnCode should be(ContinueOnReturnCodeFlag(true))

    val echoWithFalseContinueOnReturnCodeIndex =
      namespaceWithContinueOnReturnCode.workflow.calls indexWhere { call =>
        call.unqualifiedName == "echoWithFalseContinueOnReturnCode"
      }
    echoWithFalseContinueOnReturnCodeIndex should be >= 0
    namespaceWithContinueOnReturnCode.toCromwellRuntimeAttribute(echoWithFalseContinueOnReturnCodeIndex).continueOnReturnCode should be(ContinueOnReturnCodeFlag(false))

    val echoWithTrueStringContinueOnReturnCodeIndex =
      namespaceWithContinueOnReturnCode.workflow.calls indexWhere { call =>
        call.unqualifiedName == "echoWithTrueStringContinueOnReturnCode"
      }
    echoWithTrueStringContinueOnReturnCodeIndex should be >= 0
    namespaceWithContinueOnReturnCode.toCromwellRuntimeAttribute(echoWithTrueStringContinueOnReturnCodeIndex).continueOnReturnCode should be(ContinueOnReturnCodeFlag(true))

    val echoWithFalseStringContinueOnReturnCodeIndex =
      namespaceWithContinueOnReturnCode.workflow.calls indexWhere { call =>
        call.unqualifiedName == "echoWithFalseStringContinueOnReturnCode"
      }
    echoWithFalseStringContinueOnReturnCodeIndex should be >= 0
    namespaceWithContinueOnReturnCode.toCromwellRuntimeAttribute(echoWithFalseStringContinueOnReturnCodeIndex).continueOnReturnCode should be(ContinueOnReturnCodeFlag(false))
  }

  "WDL file with Googly config" should "parse up properly" in {
    val namespaceWithGooglyConfig = NamespaceWithWorkflow.load(WorkflowWithFullGooglyConfig)
    val calls = namespaceWithGooglyConfig.workflow.calls
    val callIndex = calls.indexWhere(call => call.unqualifiedName == "googly_task")
    callIndex should be >= 0

    val googlyCall = calls(callIndex)
    val attributes = CromwellRuntimeAttributes(googlyCall.task.runtimeAttributes, BackendType.LOCAL)
    attributes.cpu shouldBe 3
    val firstDisk = new Disk().setName("Disk1").setSizeGb(3L).setType("PERSISTENT_SSD").setAutoDelete(true)
    val secondDisk = new Disk().setName("Disk2").setSizeGb(500L).setType("PERSISTENT_HDD").setAutoDelete(true)

    val expectedDisks = Vector(firstDisk, secondDisk, CromwellRuntimeAttributes.LocalizationDisk)
    attributes.defaultDisks should contain theSameElementsAs expectedDisks

    val expectedZones = Vector("US_Metro", "US_Backwater")
    attributes.defaultZones foreach { z => expectedZones should contain (z) }

    attributes.memoryGB shouldBe 4
  }

  "WDL file with no Googly config" should "also parse up properly to defaults" in {
    val NamespaceWithoutGooglyConfig = NamespaceWithWorkflow.load(WorkflowWithoutGooglyConfig)
    val calls = NamespaceWithoutGooglyConfig.workflow.calls
    val callIndex = calls.indexWhere(call => call.unqualifiedName == "googly_task")
    callIndex should be >= 0

    val googlyCall = calls(callIndex)
    val attributes = CromwellRuntimeAttributes(googlyCall.task.runtimeAttributes, BackendType.LOCAL)
    attributes.cpu shouldBe CromwellRuntimeAttributes.Defaults.Cpu
    attributes.defaultDisks foreach { d => CromwellRuntimeAttributes.Defaults.Disk should contain (d) }
    attributes.defaultZones foreach { z => CromwellRuntimeAttributes.Defaults.Zones should contain (z) }
    attributes.memoryGB shouldBe CromwellRuntimeAttributes.Defaults.Memory
  }

  "WDL file with local disk Googly config" should "parse up properly" in {
    val NamespaceWithoutGooglyConfig = NamespaceWithWorkflow.load(WorkflowWithLocalDiskGooglyConfig)
    val calls = NamespaceWithoutGooglyConfig.workflow.calls
    val callIndex = calls.indexWhere(call => call.unqualifiedName == "googly_task")
    callIndex should be >= 0

    val googlyCall = calls(callIndex)
    val attributes = CromwellRuntimeAttributes(googlyCall.task.runtimeAttributes, BackendType.LOCAL)
    attributes.cpu shouldBe CromwellRuntimeAttributes.Defaults.Cpu

    val localHddDisk = new Disk().setName("local-disk").setSizeGb(123L).setType("PERSISTENT_HDD").setAutoDelete(true)
    attributes.defaultDisks should contain theSameElementsAs Vector(localHddDisk)
    attributes.defaultZones should contain theSameElementsAs CromwellRuntimeAttributes.Defaults.Zones
    attributes.memoryGB shouldBe CromwellRuntimeAttributes.Defaults.Memory
  }

  "WDL file with Googly config" should "issue warnings on the local backend" in {
    val namespace = NamespaceWithWorkflow.load(WorkflowWithFullGooglyConfig)
    val attributeMap = AttributeMap(namespace.workflow.calls.head.task.runtimeAttributes.attrs)
    val expectedString = "Found unsupported keys for backend 'LOCAL': cpu, defaultDisks, defaultZones, memory"
    attributeMap.unsupportedKeys(BackendType.LOCAL).head shouldBe expectedString
  }

  "WDL file without runtime section" should "not be accepted on JES backend as it has no docker" in {
    val ex = intercept[IllegalArgumentException] {
      NamespaceWithWorkflow.load(WorkflowWithoutRuntime).toCromwellRuntimeAttributes(BackendType.JES)
    }
    ex.getMessage should include ("Missing required keys in runtime configuration for backend 'JES': docker")
  }

  "WDL file with runtime section but no docker" should "not be accepted on JES backend" in {
    val ex = intercept[IllegalArgumentException] {
      val workflow = NamespaceWithWorkflow.load(WorkflowWithFailOnStderr).toCromwellRuntimeAttribute(0, BackendType.JES)
    }
    ex.getMessage should include ("Missing required keys in runtime configuration for backend 'JES': docker")
  }


  "WDL file with a seriously screwed up memory runtime" should "not parse" in {
    val ex = intercept[IllegalArgumentException] {
      val namespaceWithBorkedMemory = NamespaceWithWorkflow.load(WorkflowWithMessedUpMemory).toCromwellRuntimeAttributes(BackendType.LOCAL)
    }

    ex.getMessage should include ("should be of the form X Unit")
  }

  "WDL file with an invalid memory unit" should "say so" in {
    val ex = intercept[IllegalArgumentException] {
      val namespaceWithBorkedMemory = NamespaceWithWorkflow.load(WorkflowWithMessedUpMemoryUnit).toCromwellRuntimeAttributes(BackendType.LOCAL)
    }

    ex.getMessage should include ("is an invalid memory unit")
  }

  "WDL file with an invalid local disk" should "say so" in {
    val ex = intercept[IllegalArgumentException] {
      NamespaceWithWorkflow.load(WorkflowWithMessedUpLocalDisk).toCromwellRuntimeAttributes(BackendType.JES)
    }

    ex.getMessage should include(
      "'Disk1 123 LOCAL' should be in form 'NAME SIZE TYPE', with SIZE blank for LOCAL, otherwise SIZE in GB")
  }

  "WDL file with an invalid disk size" should "say so" in {
    val ex = intercept[IllegalArgumentException] {
      NamespaceWithWorkflow.load(WorkflowWithMessedUpDiskSize).toCromwellRuntimeAttributes(BackendType.JES)
    }

    ex.getMessage should include("123.0 not convertible to a Long")
  }

  "WDL file with an invalid disk type" should "say so" in {
    val ex = intercept[IllegalArgumentException] {
      NamespaceWithWorkflow.load(WorkflowWithMessedUpDiskType).toCromwellRuntimeAttributes(BackendType.JES)
    }

    ex.getMessage should include("Disk TYPE SDD should be one of LOCAL, SSD, HDD")
  }
}
