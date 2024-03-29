package cromwell.backend.google.pipelines.v2alpha1

import java.time.OffsetDateTime
import com.google.api.services.genomics.v2alpha1.model.{Accelerator, Disk, Event, Mount}
import cromwell.backend.google.pipelines.common.io.{DiskType, PipelinesApiAttachedDisk, PipelinesApiReferenceFilesDisk}
import cromwell.backend.google.pipelines.common.{GpuResource, MachineConstraints, PipelinesApiRuntimeAttributes}
import cromwell.core.ExecutionEvent
import cromwell.core.logging.JobLogger
import mouse.all._

import scala.language.postfixOps
import scala.util.Try

trait PipelinesUtilityConversions {
  def toAccelerator(gpuResource: GpuResource): Accelerator =
    new Accelerator().setCount(gpuResource.gpuCount.value.toLong).setType(gpuResource.gpuType.toString)
  def toMachineType(jobLogger: JobLogger)(attributes: PipelinesApiRuntimeAttributes): String =
    MachineConstraints.machineType(
      memory = attributes.memory,
      cpu = attributes.cpu,
      cpuPlatformOption = attributes.cpuPlatform,
      googleLegacyMachineSelection = attributes.googleLegacyMachineSelection,
      jobLogger = jobLogger
    )
  def toMounts(disks: Seq[PipelinesApiAttachedDisk]): List[Mount] = disks.map(toMount).toList
  def toDisks(disks: Seq[PipelinesApiAttachedDisk]): List[Disk] = disks.map(toDisk).toList
  def toMount(disk: PipelinesApiAttachedDisk): Mount = {
    val mount = new Mount()
      .setDisk(disk.name)
      .setPath(disk.mountPoint.pathAsString)
    disk match {
      case _: PipelinesApiReferenceFilesDisk =>
        mount.setReadOnly(true)
      case _ =>
        mount
    }
  }

  def toDisk(disk: PipelinesApiAttachedDisk): Disk = {
    val googleDisk = new Disk()
      .setName(disk.name)
      .setType(disk.diskType |> toV2DiskType)
      .setSizeGb(disk.sizeGb)
    disk match {
      case refDisk: PipelinesApiReferenceFilesDisk =>
        googleDisk.setSourceImage(refDisk.image)
      case _ =>
        googleDisk
    }
  }

  def toExecutionEvent(actionIndexToEventType: Map[Int, String])(event: Event): ExecutionEvent = {
    val groupingFromAction = for {
      rawValue <- Option(event.getDetails.get("actionId"))
      integerValue <- Try(Integer.valueOf(rawValue.toString)).toOption
      group <- actionIndexToEventType.get(integerValue)
    } yield group

    // There are both "Started pulling" and "Stopped pulling" events but these are confusing for metadata, especially on the
    // timing diagram. Create a single "Pulling <docker image>" grouping to absorb these events.
    def groupingFromPull: Option[String] = List("Started", "Stopped") flatMap { k =>
      Option(event.getDescription) collect {
        case d if d.startsWith(s"$k pulling") => "Pulling" + d.substring(s"$k pulling".length)
      }
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
