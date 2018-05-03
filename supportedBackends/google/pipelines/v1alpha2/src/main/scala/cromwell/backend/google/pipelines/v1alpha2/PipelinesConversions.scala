package cromwell.backend.google.pipelines.v1alpha2

import com.google.api.services.genomics.model.{Disk, LocalCopy, PipelineParameter}
import cromwell.backend.google.pipelines.common.io.PipelinesApiAttachedDisk
import cromwell.backend.google.pipelines.common.{PipelinesApiFileInput, PipelinesApiFileOutput, PipelinesApiInput, PipelinesApiLiteralInput}
import simulacrum.typeclass
import scala.language.implicitConversions

@typeclass trait ToParameter[A] {
  def toGooglePipelineParameter(a: A): PipelineParameter
  def toGoogleRunParameter(a: A): String
}

object PipelinesConversions {
  implicit val fileInputTo: ToParameter[PipelinesApiFileInput] = new ToParameter[PipelinesApiFileInput] {
    override def toGooglePipelineParameter(input: PipelinesApiFileInput): PipelineParameter = {
      new PipelineParameter().setName(input.name).setLocalCopy(
        new LocalCopy().setDisk(input.mount.name).setPath(input.local.pathAsString)
      )
    }

    override def toGoogleRunParameter(input: PipelinesApiFileInput): String = input.cloudPath
  }

  implicit val literalInputTo: ToParameter[PipelinesApiLiteralInput] = new ToParameter[PipelinesApiLiteralInput] {
    override def toGooglePipelineParameter(input: PipelinesApiLiteralInput): PipelineParameter = new PipelineParameter().setName(input.name)
    override def toGoogleRunParameter(input: PipelinesApiLiteralInput): String = input.value
  }

  implicit val fileOutputTo: ToParameter[PipelinesApiFileOutput] = new ToParameter[PipelinesApiFileOutput] {
    override def toGooglePipelineParameter(output: PipelinesApiFileOutput): PipelineParameter = {
      new PipelineParameter().setName(output.name).setLocalCopy(
        new LocalCopy().setDisk(output.mount.name).setPath(output.local.pathAsString)
      )
    }

    override def toGoogleRunParameter(output: PipelinesApiFileOutput): String = output.cloudPath
  }

  implicit val inputTo: ToParameter[PipelinesApiInput] = new ToParameter[PipelinesApiInput] {
    override def toGooglePipelineParameter(input: PipelinesApiInput): PipelineParameter = input match {
      case file: PipelinesApiFileInput => fileInputTo.toGooglePipelineParameter(file)
      case literal: PipelinesApiLiteralInput => literalInputTo.toGooglePipelineParameter(literal)
    }

    override def toGoogleRunParameter(input: PipelinesApiInput): String = input match {
      case file: PipelinesApiFileInput => fileInputTo.toGoogleRunParameter(file)
      case literal: PipelinesApiLiteralInput => literalInputTo.toGoogleRunParameter(literal)
    }
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
