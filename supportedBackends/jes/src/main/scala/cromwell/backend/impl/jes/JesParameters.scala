package cromwell.backend.impl.jes

import java.nio.file.Path

import com.google.api.services.genomics.model.{LocalCopy, PipelineParameter}
import cromwell.backend.impl.jes.io.JesAttachedDisk

sealed trait JesParameter {
  def name: String
  def toGooglePipelineParameter: PipelineParameter
  def toGoogleRunParameter: String
}

sealed trait JesInput extends JesParameter

final case class JesFileInput(name: String, gcs: String, local: Path, mount: JesAttachedDisk) extends JesInput {
  def toGooglePipelineParameter = {
    new PipelineParameter().setName(name).setLocalCopy(
      new LocalCopy().setDisk(mount.name).setPath(local.toString)
    )
  }
  val toGoogleRunParameter: String = gcs
  def containerPath: Path = mount.mountPoint.resolve(local)
}

final case class JesLiteralInput(name: String, value: String) extends JesInput {
  def toGooglePipelineParameter = new PipelineParameter().setName(name)
  val toGoogleRunParameter: String = value
}

final case class JesFileOutput(name: String, gcs: String, local: Path, mount: JesAttachedDisk) extends JesParameter {
  def toGooglePipelineParameter = {
    new PipelineParameter().setName(name).setLocalCopy(
      new LocalCopy().setDisk(mount.name).setPath(local.toString)
    )
  }
  val toGoogleRunParameter: String = gcs
}
