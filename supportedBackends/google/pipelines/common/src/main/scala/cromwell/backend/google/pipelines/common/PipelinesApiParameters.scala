package cromwell.backend.google.pipelines.common

import cromwell.backend.google.pipelines.common.io.JesAttachedDisk
import cromwell.core.path.Path

sealed trait JesParameter {
  def name: String
  def toGoogleRunParameter: String
}

sealed trait JesInput extends JesParameter

final case class JesFileInput(name: String, gcs: String, local: Path, mount: JesAttachedDisk) extends JesInput {
  val toGoogleRunParameter: String = gcs

  def containerPath: Path = mount.mountPoint.resolve(local)
}

final case class JesLiteralInput(name: String, value: String) extends JesInput {
  val toGoogleRunParameter: String = value
}

final case class JesFileOutput(name: String, gcs: String, local: Path, mount: JesAttachedDisk) extends JesParameter {

  val toGoogleRunParameter: String = gcs
}
