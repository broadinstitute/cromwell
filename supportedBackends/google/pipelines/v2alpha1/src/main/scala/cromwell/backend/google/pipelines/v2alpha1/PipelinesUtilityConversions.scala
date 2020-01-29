package cromwell.backend.google.pipelines.v2alpha1

import java.time.OffsetDateTime

import com.google.api.services.genomics.v2alpha1.model.{Accelerator, Disk, Event, Mount}
import cromwell.backend.google.pipelines.common.api.PipelinesApiRequestFactory.CreatePipelineParameters
import cromwell.backend.google.pipelines.common.io.{DiskType, PipelinesApiAttachedDisk}
import cromwell.backend.google.pipelines.common.{GpuResource, PipelinesApiRuntimeAttributes}
import cromwell.core.ExecutionEvent
import cromwell.core.logging.JobLogger
import mouse.all._

import scala.language.postfixOps
import scala.util.Try

trait PipelinesUtilityConversions {
  def toAccelerator(gpuResource: GpuResource) = new Accelerator().setCount(gpuResource.gpuCount.value.toLong).setType(gpuResource.gpuType.toString)
  def toMachineType(jobLogger: JobLogger)(attributes: PipelinesApiRuntimeAttributes) = MachineConstraints.machineType(attributes.memory, attributes.cpu, jobLogger)
  def toMounts(parameters: CreatePipelineParameters): List[Mount] = parameters.adjustedSizeDisks.map(toMount).toList
  def toDisks(parameters: CreatePipelineParameters): List[Disk] = parameters.adjustedSizeDisks.map(toDisk).toList
  def toMount(disk: PipelinesApiAttachedDisk) = new Mount()
    .setDisk(disk.name)
    .setPath(disk.mountPoint.pathAsString)
  def toDisk(disk: PipelinesApiAttachedDisk) = new Disk()
    .setName(disk.name)
    .setSizeGb(disk.sizeGb)
    .setType(disk.diskType |> toV2DiskType)

  def toExecutionEvent(actionIndexToEventType: Map[Int, String])(event: Event): ExecutionEvent = {
    val groupingFromAction = for {
      rawValue <- Option(event.getDetails.get("actionId"))
      integerValue <- Try(Integer.valueOf(rawValue.toString)).toOption
      group <- actionIndexToEventType.get(integerValue)
    } yield group

    // There are both "Started pulling" and "Stopped pulling" events but these are confusing for metadata, especially on the
    // timing diagram. Create a single "Pulling <docker image>" grouping to absorb these events.
    def groupingFromPull: Option[String] = List("Started", "Stopped") flatMap { k =>
      Option(event.getDescription) collect { case d if d.startsWith(s"$k pulling") => "Pulling" + d.substring(s"$k pulling".length)}
    } headOption

    ExecutionEvent(
      name = event.getDescription,
      offsetDateTime = OffsetDateTime.parse(event.getTimestamp),
      grouping = groupingFromAction.orElse(groupingFromPull)
    )
  }

  private def toV2DiskType(diskType: DiskType) = diskType match {
    case DiskType.HDD => "pd-standard"
    case DiskType.SSD => "pd-ssd"
    case DiskType.LOCAL => "local-ssd"
  }
}
