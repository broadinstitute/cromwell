package cromwell.services.metadata

import cats.syntax.validated._
import common.validation.ErrorOr.ErrorOr
import MetadataArchiveStatus._

sealed trait MetadataArchiveStatus {
  final def isArchived = this match {
    case Archived | ArchivedAndPurged => true
    case ArchiveFailed | Unarchived => false
  }
}

object MetadataArchiveStatus {

  lazy val MetadataArchiveStatusValues = Seq(Unarchived, Archived, ArchivedAndPurged, ArchiveFailed)

  def toDatabaseValue(status: MetadataArchiveStatus): Option[String] = status match {
    case Unarchived => None
    case other => Option(other.toString)
  }

  def fromDatabaseValue(status: Option[String]): ErrorOr[MetadataArchiveStatus] = status match {
    case None => Unarchived.validNel
    case Some(other) => withName(other)
  }

  def withName(str: String): ErrorOr[MetadataArchiveStatus] = MetadataArchiveStatusValues.find(_.toString.equalsIgnoreCase(str)) match {
    case Some(value) => value.validNel
    case None => s"No such MetadataArchiveStatus: $str".invalidNel
  }

  case object Unarchived extends MetadataArchiveStatus
  case object Archived extends MetadataArchiveStatus
  case object ArchivedAndPurged extends MetadataArchiveStatus // `purged` means that original data is deleted from METADATA_ENTRY table
  case object ArchiveFailed extends MetadataArchiveStatus

}
