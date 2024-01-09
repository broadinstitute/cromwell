package cromwell.backend.google.batch.models

import akka.http.scaladsl.model.ContentType
import cromwell.backend.google.batch.io.GcpBatchAttachedDisk
import cromwell.core.path.Path

import scala.concurrent.duration.FiniteDuration

sealed trait BatchParameter {
  def name: String

  def mount: GcpBatchAttachedDisk

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

  /**
    * True if this parameter represents a file; false if it represents a directory.
    */
  def isFileParameter: Boolean = this match {
    case _: GcpBatchFileInput => true
    case _: GcpBatchFileOutput => true
    case _: GcpBatchDirectoryInput => false
    case _: GcpBatchDirectoryOutput => false
  }
}

sealed trait GcpBatchInput extends BatchParameter
sealed trait GcpBatchOutput extends BatchParameter {
  def contentType: Option[ContentType] = None
}

final case class GcpBatchFileInput(name: String, cloudPath: Path, relativeHostPath: Path, mount: GcpBatchAttachedDisk)
    extends GcpBatchInput

final case class GcpBatchDirectoryInput(name: String,
                                        cloudPath: Path,
                                        relativeHostPath: Path,
                                        mount: GcpBatchAttachedDisk
) extends GcpBatchInput

final case class GcpBatchFileOutput(name: String,
                                    cloudPath: Path,
                                    relativeHostPath: Path,
                                    mount: GcpBatchAttachedDisk,
                                    optional: Boolean,
                                    secondary: Boolean,
                                    uploadPeriod: Option[FiniteDuration] = None,
                                    override val contentType: Option[ContentType] = None
) extends GcpBatchOutput

final case class GcpBatchDirectoryOutput(name: String,
                                         cloudPath: Path,
                                         relativeHostPath: Path,
                                         mount: GcpBatchAttachedDisk,
                                         optional: Boolean,
                                         secondary: Boolean,
                                         override val contentType: Option[ContentType] = None
) extends GcpBatchOutput

// TODO: Remove when support for V1 is stopped, this is only used to pass the extra_param auth file
final case class GcpBatchLiteralInput(name: String, value: String)
