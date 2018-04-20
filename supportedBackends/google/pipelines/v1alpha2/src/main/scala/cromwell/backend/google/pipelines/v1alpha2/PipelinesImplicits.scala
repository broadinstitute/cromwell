package cromwell.backend.google.pipelines.v1alpha2

import com.google.api.services.genomics.model.{Disk, LocalCopy, PipelineParameter}
import cromwell.backend.google.pipelines.common.io.JesAttachedDisk
import cromwell.backend.google.pipelines.common.{JesFileInput, JesFileOutput, JesInput, JesLiteralInput}

object PipelinesImplicits {
  implicit class EnhancedInput(val input: JesInput) extends AnyVal {
    def toGooglePipelineParameter = input match {
      case file: JesFileInput => file.toGooglePipelineParameter
      case literal: JesLiteralInput => literal.toGooglePipelineParameter
    }
  }

  implicit class EnhancedFileInput(val input: JesFileInput) extends AnyVal {
    def toGooglePipelineParameter = {
      new PipelineParameter().setName(input.name).setLocalCopy(
        new LocalCopy().setDisk(input.mount.name).setPath(input.local.pathAsString)
      )
    }
  }

  implicit class EnhancedLiteralInput(val input: JesLiteralInput) extends AnyVal {
    def toGooglePipelineParameter = new PipelineParameter().setName(input.name)
  }

  implicit class EnhancedOutput(val output: JesFileOutput) extends AnyVal {
    def toGooglePipelineParameter = {
      new PipelineParameter().setName(output.name).setLocalCopy(
        new LocalCopy().setDisk(output.mount.name).setPath(output.local.pathAsString)
      )
    }
  }

  implicit class EnhancedDisk(val disk: JesAttachedDisk) extends AnyVal {
    def toGoogleDisk: Disk = {
      new Disk().setName(disk.name)
        .setType(disk.diskType.googleTypeName)
        .setAutoDelete(true)
        .setSizeGb(disk.sizeGb)
        .setMountPoint(disk.mountPoint.toAbsolutePath.pathAsString)
    }
  }
}
