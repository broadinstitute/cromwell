package cromwell.engine.backend

import java.nio.file.Paths

import cromwell.CromwellTestkitSpec
import cromwell.CromwellTestkitSpec.TestWorkflowManagerSystem
import cromwell.backend.impl.jes.io.{DiskType, JesAttachedDisk, JesEmptyMountedDisk, JesWorkingDisk}
import cromwell.backend.validation.{ContinueOnReturnCodeFlag, ContinueOnReturnCodeSet}
import cromwell.core.WorkflowContext
import cromwell.engine.backend.jes.JesBackend
import cromwell.engine.backend.local.LocalBackend
import cromwell.engine.backend.runtimeattributes.CromwellRuntimeAttributes
import cromwell.engine.workflow.BackendCallKey
import cromwell.util.SampleWdl
import org.scalatest.prop.TableDrivenPropertyChecks._
import org.scalatest.prop.Tables.Table
import org.scalatest.{EitherValues, FlatSpec, Matchers}
import wdl4s.types.{WdlArrayType, WdlIntegerType}
import wdl4s.values.{WdlArray, WdlBoolean, WdlInteger, WdlString}

class CromwellRuntimeAttributeSpec extends FlatSpec with Matchers with EitherValues with WorkflowDescriptorBuilder {
  val workflowManagerSystem = new TestWorkflowManagerSystem
  override implicit val actorSystem = workflowManagerSystem.actorSystem
  val localBackend = new LocalBackend(CromwellTestkitSpec.DefaultLocalBackendConfigEntry, workflowManagerSystem.actorSystem)
  val jesBackend = new JesBackend(CromwellTestkitSpec.JesBackendConfigEntry, workflowManagerSystem.actorSystem)

  private def runtimeAttributes(wdl: SampleWdl, callName: String, backend: Backend, workflowOptionsJson: String = "{}"): CromwellRuntimeAttributes = {
    val root = backend match {
      case x: JesBackend => "gs://foobar"
      case x: LocalBackend => "/wf-root"
    }

    val descriptor: WorkflowDescriptor = materializeWorkflowDescriptorFromSources(workflowSources =
      wdl.asWorkflowSources(workflowOptions = workflowOptionsJson)).copy(wfContext = new WorkflowContext(root))

    val call = descriptor.namespace.workflow.callByName(callName).get
    val coercedInputs = descriptor.namespace.coerceRawInputs(wdl.rawInputs).get
    val inputs = coercedInputs collect { case (k, v) if s"${call.fullyQualifiedName}\\.[a-zA-Z0-9_-]+".r.findFirstMatchIn(k).isDefined =>
      k.replace(s"${call.fullyQualifiedName}.", "") -> v
    }

    BackendCallJobDescriptor(descriptor.copy(backend = backend), BackendCallKey(call, None, 1), inputs).callRuntimeAttributes
  }

  it should "have reasonable defaults" in {
    val defaults = CromwellRuntimeAttributes.defaults
    defaults.docker shouldEqual None
    defaults.memoryGB shouldEqual 2
    defaults.zones shouldEqual Seq("us-central1-a")
    defaults.disks shouldEqual Seq(JesAttachedDisk.parse("local-disk 10 SSD").get)
    defaults.cpu shouldEqual 1
    defaults.continueOnReturnCode shouldEqual ContinueOnReturnCodeSet(Set(0))
    defaults.failOnStderr shouldEqual false
    defaults.preemptible shouldEqual 0
  }

  it should "properly return the 'docker' runtime attribute" in {
    runtimeAttributes(SampleWdl.WorkflowWithStaticRuntime, "cgrep", localBackend).docker shouldEqual Some("ubuntu:latest")
    runtimeAttributes(SampleWdl.WorkflowWithStaticRuntime, "ps", localBackend).docker shouldEqual None
  }

  it should "reject a task on the JES backend without a docker container specified" in {
    val ex = intercept[IllegalArgumentException] {
      runtimeAttributes(SampleWdl.WorkflowWithoutRuntime, "hello", jesBackend)
    }
    ex.getMessage should include ("Missing required keys in runtime configuration for backend 'JES': docker")
  }

  it should "reject a task on the JES backend without a docker container specified (2)" in {
    val ex = intercept[IllegalArgumentException] {
      runtimeAttributes(SampleWdl.WorkflowWithoutRuntime, "hello", jesBackend)
    }
    ex.getMessage should include ("Missing required keys in runtime configuration for backend 'JES': docker")
  }

  it should "reject a task on the JES backend without a docker container specified (3)" in {
    val ex = intercept[IllegalArgumentException] {
      runtimeAttributes(SampleWdl.WorkflowWithFailOnStderr, "echoWithFailOnStderr", jesBackend)
    }
    ex.getMessage should include ("Missing required keys in runtime configuration for backend 'JES': docker")
  }

  it should "properly return the 'failOnStderr' runtime attribute" in {
    val echoWithFailOnStderr = runtimeAttributes(SampleWdl.WorkflowWithFailOnStderr, "echoWithFailOnStderr", localBackend)
    val echoWithoutFailOnStderr = runtimeAttributes(SampleWdl.WorkflowWithFailOnStderr, "echoWithoutFailOnStderr", localBackend)

    echoWithFailOnStderr.failOnStderr shouldEqual true
    echoWithoutFailOnStderr.failOnStderr shouldEqual false
  }

  it should "properly return the 'continueOnReturnCode' runtime attribute" in {
    val table = Table(
      ("callName", "continueOnRcValue"),
      ("echoWithSingleContinueOnReturnCode", ContinueOnReturnCodeSet(Set(123))),
      ("echoWithExpressionContinueOnReturnCode", ContinueOnReturnCodeSet(Set(444))),
      ("echoWithListContinueOnReturnCode", ContinueOnReturnCodeSet(Set(0,1,2,3))),
      ("echoWithTrueContinueOnReturnCode", ContinueOnReturnCodeFlag(true)),
      ("echoWithFalseContinueOnReturnCode", ContinueOnReturnCodeFlag(false)),
      ("echoWithTrueStringContinueOnReturnCode", ContinueOnReturnCodeFlag(true)),
      ("echoWithFalseStringContinueOnReturnCode", ContinueOnReturnCodeFlag(false))
    )

    forAll(table) { (callName, continueOnRcValue) =>
      runtimeAttributes(SampleWdl.WorkflowWithContinueOnReturnCode, callName, localBackend).continueOnReturnCode should be (continueOnRcValue)
    }
  }

  it should "properly return the 'cpu', 'disks', 'zones', and 'memory' attributes for a task run on JES" in {
    val attributes = runtimeAttributes(SampleWdl.WorkflowWithFullGooglyConfig, "googly_task", jesBackend)
    attributes.cpu shouldBe 3
    attributes.disks shouldEqual Vector(
      JesEmptyMountedDisk(DiskType.SSD, 3, Paths.get("/mnt/some-dir")),
      JesEmptyMountedDisk(DiskType.HDD, 500, Paths.get("/tmp/mnt")),
      CromwellRuntimeAttributes.DefaultJesWorkingDisk
    )
    attributes.zones shouldEqual Vector("US_Metro", "US_Backwater")
    attributes.memoryGB shouldBe 4
  }

  it should "properly return the 'cpu', 'disks', 'zones', and 'memory' attributes for a task run on JES (2)" in {
    val attributes = runtimeAttributes(SampleWdl.WorkflowWithLocalDiskGooglyConfig, "googly_task", jesBackend)
    val defaults = CromwellRuntimeAttributes.defaults
    attributes.cpu shouldBe defaults.cpu
    attributes.disks shouldEqual Vector(
      JesWorkingDisk(DiskType.HDD, 123)
    )
    attributes.zones shouldEqual defaults.zones
    attributes.memoryGB shouldBe defaults.memoryGB
  }

  it should "fallback to defaults" in {
    val attributes = runtimeAttributes(SampleWdl.WorkflowWithoutGooglyConfig, "googly_task", localBackend)
    val defaults = CromwellRuntimeAttributes.defaults
    attributes.cpu shouldBe defaults.cpu
    attributes.disks shouldEqual defaults.disks
    attributes.zones shouldEqual defaults.zones
    attributes.memoryGB shouldBe defaults.memoryGB
  }

  it should "detect unsupported runtime attributes on the local backend" in {
    val descriptor = materializeWorkflowDescriptorFromSources(workflowSources = SampleWdl.WorkflowWithFullGooglyConfig.asWorkflowSources())
    val call = descriptor.namespace.workflow.callByName("googly_task").get
    val unsupported = CromwellRuntimeAttributes.unsupportedKeys(call.task.runtimeAttributes.attrs.keys, BackendType.LOCAL)
    unsupported shouldEqual Set("disks", "cpu", "zones", "memory")
  }

  it should "reject a task with an invalid 'memory' attribute" in {
    val ex = intercept[IllegalArgumentException] {
      runtimeAttributes(SampleWdl.WorkflowWithMessedUpMemory, "messed_up_memory", jesBackend)
    }
    ex.getMessage should include ("should be of the form 'X Unit'")
  }

  it should "reject a task with an invalid 'memory' attribute (2)" in {
    val ex = intercept[IllegalArgumentException] {
      runtimeAttributes(SampleWdl.WorkflowWithMessedUpMemoryUnit, "messed_up_memory", jesBackend)
    }
    ex.getMessage should include ("is an invalid memory unit")
  }

  it should "reject a task with an invalid 'disks' parameter" in {
    val ex = intercept[IllegalArgumentException] {
      runtimeAttributes(SampleWdl.WorkflowWithMessedUpLocalDisk, "messed_up_disk", jesBackend)
    }
    ex.getMessage should include(
      "Disk strings should be of the format 'local-disk SIZE TYPE' or '/mount/point SIZE TYPE'"
    )
  }

  it should "reject a task with an invalid 'disks' parameter (2)" in {
    val ex = intercept[IllegalArgumentException] {
      runtimeAttributes(SampleWdl.WorkflowWithMessedUpDiskSize, "messed_up_disk", jesBackend)
    }
    ex.getMessage should include(
      "Disk strings should be of the format 'local-disk SIZE TYPE' or '/mount/point SIZE TYPE'"
    )
  }

  it should "reject a task with an invalid 'disks' parameter (3)" in {
    val ex = intercept[IllegalArgumentException] {
      runtimeAttributes(SampleWdl.WorkflowWithMessedUpDiskType, "messed_up_disk", jesBackend)
    }
    ex.getMessage should include("Disk TYPE SDD should be one of LOCAL, SSD, HDD")
  }

  it should "allow runtime attributes to be expressions that reference task inputs" in {
    val attributes = runtimeAttributes(SampleWdl.WorkflowWithRuntimeAttributeExpressions, "test", jesBackend)
    attributes.memoryGB shouldEqual 7
    attributes.disks shouldEqual Vector(
      JesEmptyMountedDisk(DiskType.SSD, 9, Paths.get("/mnt")),
      CromwellRuntimeAttributes.DefaultJesWorkingDisk
    )
    attributes.docker shouldEqual Some("ubuntu:latest")
  }

  it should "allow runtime attributes to be expressions that reference task inputs (2)" in {
    val attributes = runtimeAttributes(SampleWdl.WorkflowWithRuntimeAttributeExpressions2, "x", jesBackend)
    attributes.memoryGB shouldEqual 5
  }

  it should "allow workflow options to specify defaults for all tasks in a WDL file" in {
    val workflowOptions =
      """{
        |  "defaultRuntimeOptions": {
        |    "docker": "ubuntu:latest",
        |    "memory": "5000000 KB"
        |  }
        |}
      """.stripMargin

    Seq("x", "y", "z") foreach { taskName =>
      val attributes = runtimeAttributes(SampleWdl.WorkflowWithThreeTasksAndNoRuntime, taskName, jesBackend, workflowOptions)
      attributes.memoryGB shouldEqual 5
      attributes.docker shouldEqual Some("ubuntu:latest")
    }

    Seq("x", "y", "z") foreach { taskName =>
      val attributes = runtimeAttributes(SampleWdl.WorkflowWithThreeTasksAndNoRuntime, taskName, localBackend)
      attributes.memoryGB shouldEqual 2
      attributes.docker shouldEqual None
    }
  }

  it should "NOT allow workflow options to override values specified in a task runtime section" in {
    val workflowOptions =
      """{
        |  "defaultRuntimeOptions": {
        |    "docker": "ubuntu:latest"
        |  }
        |}
      """.stripMargin

    Seq("x", "y", "z") foreach { taskName =>
      val attributes = runtimeAttributes(SampleWdl.WorkflowWithThreeTasksWithRuntimeSections, taskName, localBackend, workflowOptions)
      attributes.docker shouldEqual Some("python:2.7")
    }

    Seq("x", "y", "z") foreach { taskName =>
      val attributes = runtimeAttributes(SampleWdl.WorkflowWithThreeTasksWithRuntimeSections, taskName, localBackend)
      attributes.docker shouldEqual Some("python:2.7")
    }
  }

  it should "allow override of 'preemptible' via workflow options" in {
    val twoAttemptsWfOptions =
      """
        |{
        |  "defaultRuntimeOptions": {
        |    "preemptible": 2
        |  }
        |}
      """.stripMargin

    runtimeAttributes(SampleWdl.WorkflowWithFourPreemptibleRetries, "x", jesBackend, twoAttemptsWfOptions).preemptible shouldEqual 4
    runtimeAttributes(SampleWdl.WorkflowWithStaticRuntime, "cgrep", jesBackend, twoAttemptsWfOptions).preemptible shouldEqual 2
    runtimeAttributes(SampleWdl.WorkflowWithFourPreemptibleRetries, "x", jesBackend).preemptible shouldEqual 4
    runtimeAttributes(SampleWdl.WorkflowWithStaticRuntime, "cgrep", jesBackend).preemptible shouldEqual 0
  }

  it should "contain only supported keys in the attributes map and coerce values to supported WdlTypes" in {
    val fullWfOptions =
      """
        |{
        |  "defaultRuntimeOptions": {
        |    "continueOnReturnCode": [0, 1, 2, 3],
        |    "bootDiskSizeGb": 10,
        |    "failOnStderr": true,
        |    "zones": "us-central1-a",
        |    "disks": "local-disk 10 SSD",
        |    "memory": "5000000 KB",
        |    "preemptible": 2,
        |    "cpu": 3
        |  }
        |}
      """.stripMargin

    val expectedJesKeys = Set("docker", "continueOnReturnCode", "bootDiskSizeGb", "failOnStderr", "zones", "disks", "memory", "preemptible", "cpu")
    val expectedLocalKeys = Set("docker", "continueOnReturnCode", "failOnStderr")

    val jesRA = runtimeAttributes(SampleWdl.WorkflowWithStaticRuntime, "cgrep", jesBackend, fullWfOptions).attributes
    jesRA.keySet should contain theSameElementsAs expectedJesKeys
    jesRA("continueOnReturnCode").valueString shouldBe WdlArray(WdlArrayType(WdlIntegerType), Seq(0, 1, 2, 3).map(WdlInteger(_))).valueString
    jesRA("failOnStderr").valueString shouldBe WdlBoolean(true).valueString
    jesRA("bootDiskSizeGb").valueString shouldBe WdlInteger(10).valueString
    jesRA("zones").valueString shouldBe WdlString("us-central1-a").valueString
    jesRA("memory").valueString shouldBe WdlString("5000000 KB").valueString
    jesRA("preemptible").valueString shouldBe WdlInteger(2).valueString
    jesRA("cpu").valueString shouldBe WdlInteger(3).valueString
    jesRA("docker").valueString shouldBe WdlString("ubuntu:latest").valueString

    val localRA = runtimeAttributes(SampleWdl.WorkflowWithStaticRuntime, "cgrep", localBackend, fullWfOptions).attributes
    localRA.keySet should contain theSameElementsAs expectedLocalKeys
    localRA("docker").valueString shouldBe WdlString("ubuntu:latest").valueString
    localRA("continueOnReturnCode").valueString shouldBe WdlArray(WdlArrayType(WdlIntegerType), Seq(0, 1, 2, 3).map(WdlInteger(_))).valueString
    localRA("failOnStderr").valueString shouldBe WdlBoolean(true).valueString
  }
}
