package cromwell.backend.google.pipelines.v2alpha1

import java.time.OffsetDateTime

import com.google.api.services.genomics.v2alpha1.model.{Accelerator, Disk, Event, Mount}
import cromwell.backend.google.pipelines.common.api.PipelinesApiRequestFactory.CreatePipelineParameters
import cromwell.backend.google.pipelines.common.io.PipelinesApiAttachedDisk
import cromwell.backend.google.pipelines.common.{GpuResource, PipelinesApiRuntimeAttributes}
import cromwell.backend.google.pipelines.v2alpha1.PipelinesConversions._
import cromwell.core.ExecutionEvent
import cromwell.core.logging.JobLogger
import simulacrum.typeclass

import scala.language.implicitConversions

@typeclass trait EventConversion[A <: Event] {
  def toExecutionEvent(event: A) : ExecutionEvent = ExecutionEvent(event.getDescription, OffsetDateTime.parse(event.getTimestamp))
}

@typeclass trait DiskConversion[A <: PipelinesApiAttachedDisk] {
  def toMount(disk: A) = new Mount()
    .setDisk(disk.name)
    .setPath(disk.mountPoint.pathAsString)

  def toDisk(disk: A) = new Disk()
    .setName(disk.name)
    .setSizeGb(disk.sizeGb)
    .setType(disk.diskType.googleTypeName)
}

@typeclass trait CreatePipelineParametersConversion[A <: CreatePipelineParameters] {
  def toMounts(parameters: A): List[Mount] = parameters.runtimeAttributes.disks.map(diskConversion.toMount).toList
  def toDisks(parameters: A): List[Disk] = parameters.runtimeAttributes.disks.map(diskConversion.toDisk).toList
}

@typeclass trait AttributesConversion[A <: PipelinesApiRuntimeAttributes] {
  def toMachineType(attributes: A, jobLogger: JobLogger) = MachineConstraints.machineType(attributes.memory, attributes.cpu, jobLogger)
}

@typeclass trait GpuResourceConversion[A <: GpuResource] {
  def toAccelerator(gpuResource: A) = new Accelerator().setCount(gpuResource.gpuCount.value.toLong).setType(gpuResource.gpuType.toString)
}

trait PipelinesUtilityConversions {
  implicit val eventConversion: EventConversion[Event] = new EventConversion[Event] {}
  implicit val diskConversion: DiskConversion[PipelinesApiAttachedDisk] = new DiskConversion[PipelinesApiAttachedDisk] {}
  implicit val createParametersConversion: CreatePipelineParametersConversion[CreatePipelineParameters] = new CreatePipelineParametersConversion[CreatePipelineParameters] {}
  implicit val attributesConversion: AttributesConversion[PipelinesApiRuntimeAttributes] = new AttributesConversion[PipelinesApiRuntimeAttributes] {}
  implicit val gpuResourceConversion: GpuResourceConversion[GpuResource] = new GpuResourceConversion[GpuResource] {}
}
