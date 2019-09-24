package cromwell.services.metadata

sealed trait MetadataArchiveStatus

object MetadataArchiveStatus {

  def toDatabaseValue(status: MetadataArchiveStatus): Option[String] = status match {
    case Unarchived => None
    case other => Option(other.toString)
  }

  case object Unarchived extends MetadataArchiveStatus
  case object Archived extends MetadataArchiveStatus
  case object Failed extends MetadataArchiveStatus

}
