package cromwell.backend.google.pipelines.v2alpha1

import java.time.OffsetDateTime

import com.google.api.services.genomics.v2alpha1.model.{Accelerator, Disk, Event, Mount}
import cromwell.backend.google.pipelines.common.api.PipelinesApiRequestFactory.CreatePipelineParameters
import cromwell.backend.google.pipelines.common.io.{DiskType, PipelinesApiAttachedDisk}
import cromwell.backend.google.pipelines.common.{GpuResource, PipelinesApiRuntimeAttributes}
import cromwell.core.ExecutionEvent
import cromwell.core.logging.JobLogger
import mouse.all._

trait PipelinesUtilityConversions {
  def toAccelerator(gpuResource: GpuResource) = new Accelerator().setCount(gpuResource.gpuCount.value.toLong).setType(gpuResource.gpuType.toString)
  def toMachineType(jobLogger: JobLogger)(attributes: PipelinesApiRuntimeAttributes) = MachineConstraints.machineType(attributes.memory, attributes.cpu, jobLogger)
  def toMounts(parameters: CreatePipelineParameters): List[Mount] = parameters.runtimeAttributes.disks.map(toMount).toList
  def toDisks(parameters: CreatePipelineParameters): List[Disk] = parameters.runtimeAttributes.disks.map(toDisk).toList
  def toMount(disk: PipelinesApiAttachedDisk) = new Mount()
    .setDisk(disk.name)
    .setPath(disk.mountPoint.pathAsString)
  def toDisk(disk: PipelinesApiAttachedDisk) = new Disk()
    .setName(disk.name)
    .setSizeGb(disk.sizeGb)
    .setType(disk.diskType |> toV2DiskType)
  def toExecutionEvent(event: Event) : ExecutionEvent = ExecutionEvent(event.getDescription, OffsetDateTime.parse(event.getTimestamp))

  private def toV2DiskType(diskType: DiskType) = diskType match {
    case DiskType.HDD => "pd-standard"
    case DiskType.SSD => "pd-ssd"
    case DiskType.LOCAL => "local-ssd"
  }
}
