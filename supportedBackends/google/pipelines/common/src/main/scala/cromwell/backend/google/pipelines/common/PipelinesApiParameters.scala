package cromwell.backend.google.pipelines.common

import akka.http.scaladsl.model.ContentType
import cromwell.backend.google.pipelines.common.io.PipelinesApiAttachedDisk
import cromwell.core.path.Path

import scala.concurrent.duration.FiniteDuration

sealed trait PipelinesParameter {
  def name: String

  def mount: PipelinesApiAttachedDisk

  /**
    * The Path where the input file resides.  The backend-specific localization
    * code handles creating actions for each specific filesystem
    * implementation.
    *
    * e.g: gs://root_bucket/input_data/my_input.bam
    */
  def cloudPath: Path

  /**
    * Relative path on the host machine where the file should be localized to / delocalized from.
    * Note that the actual localization / delocalization happens in a docker container, therefore
    * [[containerPath]] should be used as the actual source / destination path when localizing / delocalizing
    *
    * e.g: root_bucket/input_data/my_input.bam
    */
  def relativeHostPath: Path

  /**
    * Path in the docker container. It must be mounted on the docker from / to its hostPath
    *
    * e.g: /cromwell_root/root_bucket/input_data/my_input.bam
    */
  def containerPath: Path = mount.mountPoint.resolve(relativeHostPath)
}

sealed trait PipelinesApiInput extends PipelinesParameter
sealed trait PipelinesApiOutput extends PipelinesParameter {
  def contentType: Option[ContentType] = None
}

final case class PipelinesApiFileInput(name: String,
                                       cloudPath: Path,
                                       relativeHostPath: Path,
                                       mount: PipelinesApiAttachedDisk) extends PipelinesApiInput

final case class PipelinesApiDirectoryInput(name: String,
                                            cloudPath: Path,
                                            relativeHostPath: Path,
                                            mount: PipelinesApiAttachedDisk) extends PipelinesApiInput

final case class PipelinesApiFileOutput(name: String,
                                        cloudPath: Path,
                                        relativeHostPath: Path,
                                        mount: PipelinesApiAttachedDisk,
                                        optional: Boolean,
                                        secondary: Boolean,
                                        uploadPeriod: Option[FiniteDuration] = None,
                                        override val contentType: Option[ContentType] = None) extends PipelinesApiOutput

final case class PipelinesApiDirectoryOutput(name: String,
                                             cloudPath: Path,
                                             relativeHostPath: Path,
                                             mount: PipelinesApiAttachedDisk,
                                             optional: Boolean,
                                             secondary: Boolean,
                                             override val contentType: Option[ContentType] = None) extends PipelinesApiOutput

// TODO: Remove when support for V1 is stopped, this is only used to pass the extra_param auth file
final case class PipelinesApiLiteralInput(name: String, value: String)
