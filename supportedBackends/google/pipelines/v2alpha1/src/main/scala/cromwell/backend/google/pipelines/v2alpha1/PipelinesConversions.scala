package cromwell.backend.google.pipelines.v2alpha1

import java.time.OffsetDateTime

import com.google.api.services.genomics.v2alpha1.model._
import cromwell.backend.google.pipelines.common.api.PipelinesApiRequestFactory.CreatePipelineParameters
import cromwell.backend.google.pipelines.common.io.PipelinesApiAttachedDisk
import cromwell.backend.google.pipelines.common._
import cromwell.backend.google.pipelines.v2alpha1.api.ActionBuilder._
import cromwell.backend.google.pipelines.v2alpha1.api.ActionFlag
import cromwell.core.ExecutionEvent
import wdl4s.parser.MemoryUnit

object PipelinesConversions {
  implicit class EnhancedEvent(val event: Event) extends AnyVal {
    def toExecutionEvent = ExecutionEvent(event.getDescription, OffsetDateTime.parse(event.getTimestamp))
  }
  
  implicit class DiskConversion(val disk: PipelinesApiAttachedDisk) extends AnyVal {
    def toMount = new Mount()
      .setDisk(disk.name)
      .setPath(disk.mountPoint.pathAsString)

    def toDisk = new Disk()
      .setName(disk.name)
      .setSizeGb(disk.sizeGb)
      .setType(disk.diskType.googleTypeName)
  }

  implicit class EnhancedCreatePipelineParameters(val parameters: CreatePipelineParameters) extends AnyVal {
    def toMounts: List[Mount] = parameters.runtimeAttributes.disks.map(_.toMount).toList
    def toDisks: List[Disk] = parameters.runtimeAttributes.disks.map(_.toDisk).toList
  }

  implicit class EnhancedFileInput(val fileInput: PipelinesApiFileInput) extends AnyVal {
    def toEnvironment = Map(fileInput.name -> fileInput.containerPath)

    def toAction(mounts: List[Mount]) = gsutil("cp", fileInput.cloudPath, fileInput.containerPath)(mounts, description = Option("localizing"))

    def toMount = {
      new Mount()
        .setDisk(fileInput.mount.name)
        .setPath(fileInput.mount.mountPoint.pathAsString)
    }
  }

  implicit class EnhancedFileOutput(val fileOutput: PipelinesApiFileOutput) extends AnyVal {
    def toEnvironment = Map(fileOutput.name -> fileOutput.containerPath)

    def toAction(mounts: List[Mount], gsutilFlags: List[String] = List.empty) = {
      gsutil("cp", fileOutput.containerPath, fileOutput.cloudPath)(mounts, List(ActionFlag.AlwaysRun), description = Option("delocalizing"))
    }

    def toMount = {
      new Mount()
        .setDisk(fileOutput.mount.name)
        .setPath(fileOutput.mount.mountPoint.pathAsString)
    }
  }

  implicit class EnhancedinputLiteral(val literalInput: PipelinesApiLiteralInput) extends AnyVal {
    def toEnvironment = Map(literalInput.name -> literalInput.value)
  }

  implicit class EnhancedAttributes(val attributes: PipelinesApiRuntimeAttributes) extends AnyVal {
    def toMachineType = {
      val cpu = attributes.cpu
      // https://cloud.google.com/genomics/reference/rpc/google.genomics.v2alpha1#virtualmachine
      // https://cloud.google.com/compute/docs/instances/creating-instance-with-custom-machine-type
      val memory = attributes.memory.to(MemoryUnit.MB).asRoundedUpMultipleOf(256).amount.toInt
      s"custom-$cpu-$memory"
    }
  }

  implicit class EnhancedGpuResource(val gpuResource: GpuResource) extends AnyVal {
    def toAccelerator = new Accelerator().setCount(gpuResource.gpuCount.toLong).setType(gpuResource.gpuType.toString)
  }
}
