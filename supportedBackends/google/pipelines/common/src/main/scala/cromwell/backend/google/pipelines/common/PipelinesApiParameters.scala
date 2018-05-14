package cromwell.backend.google.pipelines.common

import akka.http.scaladsl.model.ContentType
import cromwell.backend.google.pipelines.common.io.PipelinesApiAttachedDisk
import cromwell.core.path.Path

sealed trait PipelinesParameter {
  def name: String
}

sealed trait PipelinesApiFileParameter extends PipelinesParameter {
  /**
    * This HAS to be GCS, for now
    * The file already exists in GCS for input files, but does not for output files
    */
  def cloudPath: String

  /**
    * Path in the docker container. It must be mounted on the docker from / to its hostPath
    */
  def containerPath: Path
}

sealed trait PipelinesApiInput extends PipelinesParameter

final case class PipelinesApiFileInput(name: String,
                                       cloudPath: String,
                                       local: Path,
                                       mount: PipelinesApiAttachedDisk) extends PipelinesApiFileParameter with PipelinesApiInput {
  def containerPath: Path = mount.mountPoint.resolve(local)
}

final case class PipelinesApiFileOutput(name: String,
                                        cloudPath: String,
                                        local: Path,
                                        mount: PipelinesApiAttachedDisk,
                                        optional: Boolean,
                                        contentType: Option[ContentType] = None) extends PipelinesApiFileParameter {
  def containerPath: Path = mount.mountPoint.resolve(local)
}

final case class PipelinesApiLiteralInput(name: String, value: String) extends PipelinesApiInput
