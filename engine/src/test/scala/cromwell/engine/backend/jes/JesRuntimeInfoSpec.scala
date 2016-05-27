package cromwell.engine.backend.jes

import cromwell.engine.backend.runtimeattributes._
import org.scalatest.{FlatSpec, Matchers}

class JesRuntimeInfoSpec extends FlatSpec with Matchers {
  val disks: List[JesAttachedDisk] = List(new JesWorkingDisk(DiskType.SSD, 10))
  val zones: List[String] = List("central-1")

  it should "create Preemptible Runtime Info" in {
    val attributes = new CromwellRuntimeAttributes(attributes = Map.empty,
      docker = Some("docker"),
      zones = zones,
      failOnStderr = true,
      continueOnReturnCode = ContinueOnReturnCodeFlag(false),
      cpu = 2,
      // NOTE: This value is irrelevant for this test as it is not sufficient to determine if a Call should be started with a preemptible VM or not
      preemptible = 1,
      disks = disks,
      memoryGB = 4,
      queue = Some("workq"),
      walltime = "24:00:00",
      bootDiskSizeGb = CromwellRuntimeAttributes.defaults.bootDiskSizeGb)

    val runtimeInfo = PreemptibleJesRuntimeInfo("command", attributes)
    runtimeInfo.resources.getPreemptible shouldBe true
    runtimeInfo.resources.getMinimumCpuCores shouldBe 2
    runtimeInfo.resources.getMinimumRamGb shouldBe 4D
    runtimeInfo.resources.getDisks should contain theSameElementsAs disks.map(_.toGoogleDisk)
    runtimeInfo.resources.getZones should contain theSameElementsAs zones
    runtimeInfo.docker.getCmd shouldBe "command"
    runtimeInfo.docker.getImageName shouldBe "docker"
  }

  it should "create NonPreemptible Runtime Info" in {
    val attributes = new CromwellRuntimeAttributes(attributes = Map.empty,
      docker = Some("docker"),
      zones = zones,
      failOnStderr = true,
      continueOnReturnCode = ContinueOnReturnCodeFlag(false),
      cpu = 2,
      // NOTE: This value is irrelevant for this test as it is not sufficient to determine if a Call should be started with a preemptible VM or not
      preemptible = 3,
      disks = disks,
      memoryGB = 4,
      queue = Some("workq"),
      walltime = "24:00:00",
      bootDiskSizeGb = CromwellRuntimeAttributes.defaults.bootDiskSizeGb)

    val runtimeInfo = NonPreemptibleJesRuntimeInfo("command", attributes)
    runtimeInfo.resources.getPreemptible shouldBe false
    runtimeInfo.resources.getMinimumCpuCores shouldBe 2
    runtimeInfo.resources.getMinimumRamGb shouldBe 4D
    runtimeInfo.resources.getDisks should contain theSameElementsAs disks.map(_.toGoogleDisk)
    runtimeInfo.resources.getZones should contain theSameElementsAs zones
    runtimeInfo.docker.getCmd shouldBe "command"
    runtimeInfo.docker.getImageName shouldBe "docker"
  }

}
