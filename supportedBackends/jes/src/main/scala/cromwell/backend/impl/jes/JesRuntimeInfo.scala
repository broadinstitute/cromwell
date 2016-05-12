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
}

object NonPreemptibleJesRuntimeInfo extends JesRuntimeInfoBuilder {
  def apply(commandLine: String, runtimeAttributes: JesRuntimeAttributes): JesRuntimeInfo = {
    /*
     It should be impossible for docker to be None here. Enforcing that w/ ADTs seemed more trouble than
      it was worth. If you're ever debugging a NoSuchElementException which leads you here, that means
      the more trouble than worth calculation was incorrect and we should have separate RuntimeAttributes for
      docker and no docker cases
    */
    val dockerImage = runtimeAttributes.dockerImage.get
    val resources = buildResources(runtimeAttributes).setPreemptible(false)
    new NonPreemptibleJesRuntimeInfo(resources, buildDockerExecutor(commandLine, dockerImage))
  }
}

object PreemptibleJesRuntimeInfo extends JesRuntimeInfoBuilder {
  def apply(commandLine: String, runtimeAttributes: JesRuntimeAttributes): JesRuntimeInfo = {
    // See comment above
    val dockerImage = runtimeAttributes.dockerImage.get
    val resources = buildResources(runtimeAttributes).setPreemptible(true)
    new PreemptibleJesRuntimeInfo(resources, buildDockerExecutor(commandLine, dockerImage))
  }
}

case class NonPreemptibleJesRuntimeInfo(resources: PipelineResources, docker: DockerExecutor) extends JesRuntimeInfo
case class PreemptibleJesRuntimeInfo(resources: PipelineResources, docker: DockerExecutor) extends JesRuntimeInfo
