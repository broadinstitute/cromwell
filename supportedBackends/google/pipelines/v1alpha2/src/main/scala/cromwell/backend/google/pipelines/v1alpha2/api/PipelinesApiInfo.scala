package cromwell.backend.google.pipelines.v1alpha2.api

import com.google.api.services.genomics.model.{DockerExecutor, PipelineResources}
import cromwell.backend.google.pipelines.common.api.PipelinesApiRequestFactory.CreatePipelineParameters
import cromwell.backend.google.pipelines.common.{GpuResource, PipelinesApiRuntimeAttributes}
import cromwell.backend.google.pipelines.v1alpha2.PipelinesConversions._
import wdl4s.parser.MemoryUnit

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

  def buildResources(createPipelineParameters: CreatePipelineParameters): PipelineResources = {
    val runtimeAttributes =createPipelineParameters.runtimeAttributes
    val resources = new PipelineResources()
      .setMinimumRamGb(runtimeAttributes.memory.to(MemoryUnit.GB).amount)
      .setMinimumCpuCores(runtimeAttributes.cpu.value)
      .setZones(runtimeAttributes.zones.asJava)
      .setDisks(runtimeAttributes.disks.map(_.toGoogleDisk).asJava)
      .setBootDiskSizeGb(runtimeAttributes.bootDiskSize)
      .setNoAddress(createPipelineParameters.effectiveNoAddressValue)

      setGpu(resources, runtimeAttributes)
  }

  def build(commandLine: String, createPipelineParameters: CreatePipelineParameters): PipelineInfo
}

object NonPreemptiblePipelineInfoBuilder extends PipelineInfoBuilder {
  def build(commandLine: String, createPipelineParameters: CreatePipelineParameters): PipelineInfo = {
    val resources = buildResources(createPipelineParameters).setPreemptible(false)
    new NonPreemptiblePipelineInfoBuilder(resources, buildDockerExecutor(commandLine, createPipelineParameters.dockerImage))
  }
}

object PreemptiblePipelineInfoBuilder extends PipelineInfoBuilder {
  def build(commandLine: String, createPipelineParameters: CreatePipelineParameters): PipelineInfo = {
    val resources = buildResources(createPipelineParameters).setPreemptible(true)
    new PreemptiblePipelineInfoBuilder(resources, buildDockerExecutor(commandLine, createPipelineParameters.dockerImage))
  }
}

case class NonPreemptiblePipelineInfoBuilder private(resources: PipelineResources, docker: DockerExecutor) extends PipelineInfo
case class PreemptiblePipelineInfoBuilder private(resources: PipelineResources, docker: DockerExecutor) extends PipelineInfo
