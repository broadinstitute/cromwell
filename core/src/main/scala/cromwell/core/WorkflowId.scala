package cromwell.core

import java.util.UUID

sealed trait WorkflowId {
  def id: UUID

  override def toString = id.toString
  def shortString = id.toString.split("-")(0)

  def toRoot: RootWorkflowId =
    this match {
      case root: RootWorkflowId => root
      case _ => RootWorkflowId(id)
    }

  def toPossiblyNotRoot: PossiblyNotRootWorkflowId =
    this match {
      case possiblyNotRoot: PossiblyNotRootWorkflowId => possiblyNotRoot
      case _ => PossiblyNotRootWorkflowId(id)
    }
}

object WorkflowId {
  def apply(id: UUID): WorkflowId = PossiblyNotRootWorkflowId(id)

  def unapply(arg: WorkflowId): Option[UUID] = Option(arg.id)

  def fromString(id: String): WorkflowId = WorkflowId(UUID.fromString(id))

  def randomId(): WorkflowId = WorkflowId(UUID.randomUUID())
}

final case class RootWorkflowId(override val id: UUID) extends WorkflowId

final case class PossiblyNotRootWorkflowId(override val id: UUID) extends WorkflowId
