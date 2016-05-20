package cromwell.database.obj

import java.sql.Timestamp
import java.time.OffsetDateTime

import cromwell.core.{WorkflowMetadataKeys, WorkflowState}

import scalaz.Scalaz._
import scalaz.Semigroup
import cromwell.database.SqlConverters._

object WorkflowMetadataSummary {
  implicit val WorkflowMetadataSummarySemigroup = new Semigroup[WorkflowMetadataSummary] {
    override def append(f1: WorkflowMetadataSummary, f2: => WorkflowMetadataSummary): WorkflowMetadataSummary = f1.append(f2)
  }

  implicit class MetadatumEnhancer(val metadatum: Metadatum) extends AnyVal {
    def toSummary: WorkflowMetadataSummary = {
      val base = WorkflowMetadataSummary(metadatum.workflowUuid)
      metadatum.key match {
        case WorkflowMetadataKeys.Name => base.copy(name = metadatum.value)
        case WorkflowMetadataKeys.Status => base.copy(status = metadatum.value)
        case WorkflowMetadataKeys.StartTime => base.copy(startDate = metadatum.value map OffsetDateTime.parse map { _.toSystemTimestamp })
        case WorkflowMetadataKeys.EndTime => base.copy(endDate = metadatum.value map OffsetDateTime.parse map { _.toSystemTimestamp })
      }
    }
  }
}


case class WorkflowMetadataSummary
(
  workflowUuid: String,
  name: Option[String] = None,
  status: Option[String] = None,
  startDate: Option[Timestamp] = None,
  endDate: Option[Timestamp] = None,
  workflowMetadataSummaryId: Option[Long] = None
) {

  def append(that: WorkflowMetadataSummary): WorkflowMetadataSummary = {
    // Resolve the status if both `this` and `that` have defined statuses.  This will evaluate to `None`
    // if one or both of the statuses is not defined.
    val resolvedStatus = for {
      a <- this.status map WorkflowState.fromString
      b <- that.status map WorkflowState.fromString
    } yield (a |+| b).toString

    WorkflowMetadataSummary(
      workflowUuid = workflowUuid,
      // If not both statuses are defined, take whichever is defined.
      status = resolvedStatus orElse this.status orElse that.status,
      name = this.name orElse that.name,
      startDate = this.startDate orElse that.startDate,
      endDate = this.endDate orElse that.endDate
    )
  }
}
