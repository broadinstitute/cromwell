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
  
  def setGpu(resources: PipelineResources, runtimeAttributes: JesRuntimeAttributes) = {
    runtimeAttributes.gpuResource match {
      case Some(GpuResource(gpuType, gpuCount)) => resources
        .set("acceleratorType", gpuType)
        .set("acceleratorCount", gpuCount)
      case _ => resources
    }
  }

  def buildResources(runtimeAttributes: JesRuntimeAttributes): PipelineResources = {
    val resources = new PipelineResources()
      .setMinimumRamGb(runtimeAttributes.memory.to(MemoryUnit.GB).amount)
      .setMinimumCpuCores(runtimeAttributes.cpu)
      .setZones(runtimeAttributes.zones.asJava)
      .setDisks(runtimeAttributes.disks.map(_.toGoogleDisk).asJava)
      .setBootDiskSizeGb(runtimeAttributes.bootDiskSize)
      .setNoAddress(runtimeAttributes.noAddress)

      setGpu(resources, runtimeAttributes)
  }

  def build(commandLine: String, runtimeAttributes: JesRuntimeAttributes, docker: String): JesPipelineInfo
}

object NonPreemptibleJesPipelineInfoBuilder extends JesPipelineInfoBuilder {
  def build(commandLine: String, runtimeAttributes: JesRuntimeAttributes, dockerImage: String): JesPipelineInfo = {
    val resources = buildResources(runtimeAttributes).setPreemptible(false)
    new NonPreemptibleJesPipelineInfoBuilder(resources, buildDockerExecutor(commandLine, dockerImage))
  }
}

object PreemptibleJesPipelineInfoBuilder extends JesPipelineInfoBuilder {
  def build(commandLine: String, runtimeAttributes: JesRuntimeAttributes, dockerImage: String): JesPipelineInfo = {
    val resources = buildResources(runtimeAttributes).setPreemptible(true)
    new PreemptibleJesPipelineInfoBuilder(resources, buildDockerExecutor(commandLine, dockerImage))
  }
}

case class NonPreemptibleJesPipelineInfoBuilder private(resources: PipelineResources, docker: DockerExecutor) extends JesPipelineInfo
case class PreemptibleJesPipelineInfoBuilder private(resources: PipelineResources, docker: DockerExecutor) extends JesPipelineInfo
