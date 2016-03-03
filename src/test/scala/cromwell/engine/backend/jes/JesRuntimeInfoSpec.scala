package cromwell.engine.backend.jes

import com.google.api.services.genomics.model.Disk
import cromwell.engine.backend.runtimeattributes.{ContinueOnReturnCodeFlag, CromwellRuntimeAttributes}
import org.scalatest.{FlatSpec, Matchers}

class JesRuntimeInfoSpec extends FlatSpec with Matchers {

  it should "create Preemptible Runtime Info" in {
    val disks: List[Disk] = List(new Disk())
    val zones: List[String] = List("central-1")
    val attributes = new CromwellRuntimeAttributes(attributes = Map.empty,
      docker = Some("docker"),
      zones = zones,
      failOnStderr = true,
      continueOnReturnCode = ContinueOnReturnCodeFlag(false),
      cpu = 2L,
      // NOTE: This value is irrelevant for this test as it is not sufficient to determine if a Call should be started with a preemptible VM or not
      preemptible = 1,
      disks = disks,
      memoryGB = 4D)

    val runtimeInfo = PreemptibleJesRuntimeInfo("command", attributes)
    runtimeInfo.resources.getPreemptible shouldBe true
    runtimeInfo.resources.getCpu shouldBe 2L
    runtimeInfo.resources.getRamGb shouldBe 4D
    runtimeInfo.resources.getDisks should contain theSameElementsAs disks
    runtimeInfo.resources.getZones should contain theSameElementsAs zones
    runtimeInfo.docker.getCmd shouldBe "command"
    runtimeInfo.docker.getImage shouldBe "docker"
  }

  it should "create NonPreemptible Runtime Info" in {
    val disks: List[Disk] = List(new Disk())
    val zones: List[String] = List("central-1")
    val attributes = new CromwellRuntimeAttributes(attributes = Map.empty,
      docker = Some("docker"),
      zones = zones,
      failOnStderr = true,
      continueOnReturnCode = ContinueOnReturnCodeFlag(false),
      cpu = 2L,
      // NOTE: This value is irrelevant for this test as it is not sufficient to determine if a Call should be started with a preemptible VM or not
      preemptible = 3,
      disks = disks,
      memoryGB = 4D)

    val runtimeInfo = NonPreemptibleJesRuntimeInfo("command", attributes)
    runtimeInfo.resources.getPreemptible shouldBe false
    runtimeInfo.resources.getCpu shouldBe 2L
    runtimeInfo.resources.getRamGb shouldBe 4D
    runtimeInfo.resources.getDisks should contain theSameElementsAs disks
    runtimeInfo.resources.getZones should contain theSameElementsAs zones
    runtimeInfo.docker.getCmd shouldBe "command"
    runtimeInfo.docker.getImage shouldBe "docker"
  }

}
