package cromwell.engine.backend.jes

import com.google.api.services.genomics.model.Disk
import cromwell.engine.backend.runtimeattributes.{ContinueOnReturnCodeFlag, ContinueOnReturnCode, CromwellRuntimeAttributes}
import cromwell.util.SampleWdl.ContinueOnReturnCode
import org.scalatest.{Matchers, FlatSpec}

class JesRuntimeInfoSpec extends FlatSpec with Matchers {

  it should "create Preemptible Runtime Info" in {
    val disks: List[Disk] = List(new Disk())
    val zones: List[String] = List("central-1")
    val attributes = new CromwellRuntimeAttributes(docker = Some("docker"),
      defaultZones = zones,
      failOnStderr = true,
      continueOnReturnCode = ContinueOnReturnCodeFlag(false),
      cpu = 2L,
      // This value is irrelevant at this point as the decision of creating a pre-emptible VM or not has been done before creating the runtime info
      preemptible = None,
      defaultDisks = disks,
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
    val attributes = new CromwellRuntimeAttributes(docker = Some("docker"),
      defaultZones = zones,
      failOnStderr = true,
      continueOnReturnCode = ContinueOnReturnCodeFlag(false),
      cpu = 2L,
      // This value is irrelevant at this point as the decision of creating a pre-emptible VM or not has been done before creating the runtime info
      preemptible = Some(1),
      defaultDisks = disks,
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
