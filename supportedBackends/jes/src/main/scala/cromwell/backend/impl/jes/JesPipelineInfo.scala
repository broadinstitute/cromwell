package cromwell.backend.impl.jes

import com.google.api.services.genomics.model.{DockerExecutor, PipelineResources}
import wdl4s.parser.MemoryUnit

import scala.collection.JavaConverters._

sealed trait JesPipelineInfo {
  def resources: PipelineResources
  def docker: DockerExecutor
}

trait JesPipelineInfoBuilder {
  def buildDockerExecutor(commandLine: String, dockerImage: String): DockerExecutor = {
    val docker = new DockerExecutor()
    docker.setImageName(dockerImage).setCmd(commandLine)
  }

  def buildResources(runtimeAttributes: JesRuntimeAttributes): PipelineResources = {
    new PipelineResources()
      .setMinimumRamGb(runtimeAttributes.memory.to(MemoryUnit.GB).amount)
      .setMinimumCpuCores(runtimeAttributes.cpu)
      .setZones(runtimeAttributes.zones.asJava)
      .setDisks(runtimeAttributes.disks.map(_.toGoogleDisk).asJava)
      .setBootDiskSizeGb(runtimeAttributes.bootDiskSize)
      .set(Run.NoAddressFieldName, runtimeAttributes.noAddress)
  }

  def build(commandLine: String, runtimeAttributes: JesRuntimeAttributes): JesPipelineInfo
}

object NonPreemptibleJesPipelineInfoBuilder extends JesPipelineInfoBuilder {
  def build(commandLine: String, runtimeAttributes: JesRuntimeAttributes): JesPipelineInfo = {
    /*
     It should be impossible for docker to be None here. Enforcing that w/ ADTs seemed more trouble than
      it was worth. If you're ever debugging a NoSuchElementException which leads you here, that means
      the more trouble than worth calculation was incorrect and we should have separate RuntimeAttributes for
      docker and no docker cases
    */
    val dockerImage = runtimeAttributes.dockerImage.get
    val resources = buildResources(runtimeAttributes).setPreemptible(false)
    new NonPreemptibleJesPipelineInfoBuilder(resources, buildDockerExecutor(commandLine, dockerImage))
  }
}

object PreemptibleJesPipelineInfoBuilder extends JesPipelineInfoBuilder {
  def build(commandLine: String, runtimeAttributes: JesRuntimeAttributes): JesPipelineInfo = {
    // See comment above
    val dockerImage = runtimeAttributes.dockerImage.get
    val resources = buildResources(runtimeAttributes).setPreemptible(true)
    new PreemptibleJesPipelineInfoBuilder(resources, buildDockerExecutor(commandLine, dockerImage))
  }
}

case class NonPreemptibleJesPipelineInfoBuilder private(resources: PipelineResources, docker: DockerExecutor) extends JesPipelineInfo
case class PreemptibleJesPipelineInfoBuilder private(resources: PipelineResources, docker: DockerExecutor) extends JesPipelineInfo
