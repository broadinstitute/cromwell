package cromwell.backend.google.pipelines.v1alpha2.api

import com.google.api.services.genomics.model.{DockerExecutor, PipelineResources}
import cromwell.backend.google.pipelines.common.{GpuResource, PipelinesApiRuntimeAttributes}
import cromwell.backend.google.pipelines.v1alpha2.PipelinesConversions._

import scala.collection.JavaConverters._

sealed trait PipelineInfo {
  def resources: PipelineResources
  def docker: DockerExecutor
}

trait PipelineInfoBuilder {
  def buildDockerExecutor(commandLine: String, dockerImage: String): DockerExecutor = {
    val docker = new DockerExecutor()
    docker.setImageName(dockerImage).setCmd(commandLine)
  }
  
  def setGpu(resources: PipelineResources, runtimeAttributes: PipelinesApiRuntimeAttributes) = {
    runtimeAttributes.gpuResource match {
      case Some(GpuResource(gpuType, gpuCount, _)) => resources
        .setAcceleratorType(gpuType.toString)
        .setAcceleratorCount(gpuCount.value.toLong)
      case _ => resources
    }
  }

  def buildResources(runtimeAttributes: PipelinesApiRuntimeAttributes): PipelineResources = {
    val resources = new PipelineResources()
      .setMinimumRamGb(runtimeAttributes.memory.toGigabytes)
      .setMinimumCpuCores(runtimeAttributes.cpu.value)
      .setZones(runtimeAttributes.zones.asJava)
      .setDisks(runtimeAttributes.disks.map(_.toGoogleDisk).asJava)
      .setBootDiskSizeGb(runtimeAttributes.bootDiskSize)
      .setNoAddress(runtimeAttributes.noAddress)

      setGpu(resources, runtimeAttributes)
  }

  def build(commandLine: String, runtimeAttributes: PipelinesApiRuntimeAttributes, docker: String): PipelineInfo
}

object NonPreemptiblePipelineInfoBuilder extends PipelineInfoBuilder {
  def build(commandLine: String, runtimeAttributes: PipelinesApiRuntimeAttributes, dockerImage: String): PipelineInfo = {
    val resources = buildResources(runtimeAttributes).setPreemptible(false)
    new NonPreemptiblePipelineInfoBuilder(resources, buildDockerExecutor(commandLine, dockerImage))
  }
}

object PreemptiblePipelineInfoBuilder extends PipelineInfoBuilder {
  def build(commandLine: String, runtimeAttributes: PipelinesApiRuntimeAttributes, dockerImage: String): PipelineInfo = {
    val resources = buildResources(runtimeAttributes).setPreemptible(true)
    new PreemptiblePipelineInfoBuilder(resources, buildDockerExecutor(commandLine, dockerImage))
  }
}

case class NonPreemptiblePipelineInfoBuilder private(resources: PipelineResources, docker: DockerExecutor) extends PipelineInfo
case class PreemptiblePipelineInfoBuilder private(resources: PipelineResources, docker: DockerExecutor) extends PipelineInfo
