package cromwell.services.metadata.impl

sealed trait FileDeletionStatus

object FileDeletionStatus {

  lazy val FileDeletionStatusValues = Seq(InProgress, Succeeded, Failed)

  def toDatabaseValue(status: FileDeletionStatus): String = status.toString

  case object InProgress extends FileDeletionStatus
  case object Succeeded extends FileDeletionStatus
  case object Failed extends FileDeletionStatus
}
