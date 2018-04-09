package cromwell.backend.google.pipelines.v1alpha2

import com.google.api.services.genomics.model.{Disk, LocalCopy, PipelineParameter}
import cromwell.backend.google.pipelines.common.io.PipelinesApiAttachedDisk
import cromwell.backend.google.pipelines.common.{PipelinesApiFileInput, PipelinesApiFileOutput, PipelinesApiInput, PipelinesApiLiteralInput}

object PipelinesConversions {
  implicit class EnhancedInput(val input: PipelinesApiInput) extends AnyVal {
    def toGooglePipelineParameter: PipelineParameter = input match {
      case file: PipelinesApiFileInput => file.toGooglePipelineParameter
      case literal: PipelinesApiLiteralInput => literal.toGooglePipelineParameter
    }

    def toGoogleRunParameter: String = input match {
      case file: PipelinesApiFileInput => file.toGoogleRunParameter
      case literal: PipelinesApiLiteralInput => literal.toGoogleRunParameter
    }
  }

  implicit class EnhancedFileInput(val input: PipelinesApiFileInput) extends AnyVal {
    def toGooglePipelineParameter: PipelineParameter = {
      new PipelineParameter().setName(input.name).setLocalCopy(
        new LocalCopy().setDisk(input.mount.name).setPath(input.local.pathAsString)
      )
    }

    def toGoogleRunParameter: String = input.cloudPath
  }

  implicit class EnhancedLiteralInput(val input: PipelinesApiLiteralInput) extends AnyVal {
    def toGooglePipelineParameter: PipelineParameter = new PipelineParameter().setName(input.name)
    def toGoogleRunParameter: String = input.value
  }

  implicit class EnhancedOutput(val output: PipelinesApiFileOutput) extends AnyVal {
    def toGooglePipelineParameter: PipelineParameter = {
      new PipelineParameter().setName(output.name).setLocalCopy(
        new LocalCopy().setDisk(output.mount.name).setPath(output.local.pathAsString)
      )
    }

    def toGoogleRunParameter: String = output.cloudPath
  }

  implicit class EnhancedDisk(val disk: PipelinesApiAttachedDisk) extends AnyVal {
    def toGoogleDisk: Disk = {
      new Disk().setName(disk.name)
        .setType(disk.diskType.googleTypeName)
        .setAutoDelete(true)
        .setSizeGb(disk.sizeGb)
        .setMountPoint(disk.mountPoint.toAbsolutePath.pathAsString)
    }
  }
}
