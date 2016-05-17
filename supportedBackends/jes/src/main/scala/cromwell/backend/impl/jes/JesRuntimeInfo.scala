package cromwell.backend.impl.jes

import com.google.api.services.genomics.model.{DockerExecutor, PipelineResources}
import wdl4s.parser.MemoryUnit

import scala.collection.JavaConverters._

sealed trait JesRuntimeInfo {
  def resources: PipelineResources
  def docker: DockerExecutor
}

trait JesRuntimeInfoBuilder {
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
  }

  def build(commandLine: String, runtimeAttributes: JesRuntimeAttributes): JesRuntimeInfo
}

object NonPreemptibleJesRuntimeInfoBuilder extends JesRuntimeInfoBuilder {
  def build(commandLine: String, runtimeAttributes: JesRuntimeAttributes): JesRuntimeInfo = {
    /*
     It should be impossible for docker to be None here. Enforcing that w/ ADTs seemed more trouble than
      it was worth. If you're ever debugging a NoSuchElementException which leads you here, that means
      the more trouble than worth calculation was incorrect and we should have separate RuntimeAttributes for
      docker and no docker cases
    */
    val dockerImage = runtimeAttributes.dockerImage.get
    val resources = buildResources(runtimeAttributes).setPreemptible(false)
    new NonPreemptibleJesRuntimeInfoBuilder(resources, buildDockerExecutor(commandLine, dockerImage))
  }
}

object PreemptibleJesRuntimeInfoBuilder extends JesRuntimeInfoBuilder {
  def build(commandLine: String, runtimeAttributes: JesRuntimeAttributes): JesRuntimeInfo = {
    // See comment above
    val dockerImage = runtimeAttributes.dockerImage.get
    val resources = buildResources(runtimeAttributes).setPreemptible(true)
    new PreemptibleJesRuntimeInfoBuilder(resources, buildDockerExecutor(commandLine, dockerImage))
  }
}

case class NonPreemptibleJesRuntimeInfoBuilder private(resources: PipelineResources, docker: DockerExecutor) extends JesRuntimeInfo
case class PreemptibleJesRuntimeInfoBuilder private(resources: PipelineResources, docker: DockerExecutor) extends JesRuntimeInfo
