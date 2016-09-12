package cromwell.services.metadata.impl

import java.time.OffsetDateTime

import cromwell.core.{WorkflowId, WorkflowMetadataKeys, WorkflowState}
import cromwell.database.sql.SqlConverters._
import cromwell.database.sql.tables.{Metadatum, WorkflowMetadataSummary}
import cromwell.services.ServicesStore
import cromwell.services.metadata.MetadataService.{QueryMetadata, WorkflowQueryResponse}
import cromwell.services.metadata._

import scala.concurrent.{ExecutionContext, Future}
import scalaz.Scalaz._
import scalaz.{NonEmptyList, Semigroup}

object MetadataDatabaseAccess {

  private lazy val WorkflowMetadataSummarySemigroup = new Semigroup[WorkflowMetadataSummary] {
    override def append(summary1: WorkflowMetadataSummary,
                        summary2: => WorkflowMetadataSummary): WorkflowMetadataSummary = {
      // Resolve the status if both `this` and `that` have defined statuses.  This will evaluate to `None`
      // if one or both of the statuses is not defined.
      val resolvedStatus = for {
        thisStatus <- summary1.status map WorkflowState.fromString
        thatStatus <- summary2.status map WorkflowState.fromString
      } yield (thisStatus |+| thatStatus).toString

      WorkflowMetadataSummary(
        workflowUuid = summary1.workflowUuid,
        // If not both statuses are defined, take whichever is defined.
        status = resolvedStatus orElse summary1.status orElse summary2.status,
        name = summary1.name orElse summary2.name,
        startDate = summary1.startDate orElse summary2.startDate,
        endDate = summary1.endDate orElse summary2.endDate
      )
    }
  }

  def baseSummary(workflowUuid: String) = WorkflowMetadataSummary(workflowUuid, None, None, None, None, None)

  private implicit class MetadatumEnhancer(val metadatum: Metadatum) extends AnyVal {
    def toSummary: WorkflowMetadataSummary = {
      val base = baseSummary(metadatum.workflowUuid)
      metadatum.key match {
        case WorkflowMetadataKeys.Name => base.copy(name = metadatum.value)
        case WorkflowMetadataKeys.Status => base.copy(status = metadatum.value)
        case WorkflowMetadataKeys.StartTime =>
          base.copy(startDate = metadatum.value map OffsetDateTime.parse map { _.toSystemTimestamp })
        case WorkflowMetadataKeys.EndTime =>
          base.copy(endDate = metadatum.value map OffsetDateTime.parse map { _.toSystemTimestamp })
      }
    }
  }

  private def buildUpdatedSummary(existingSummary: Option[WorkflowMetadataSummary],
                                  metadataForUuid: Seq[Metadatum]): WorkflowMetadataSummary = {
    implicit val wmss = WorkflowMetadataSummarySemigroup

    val base = existingSummary.getOrElse(baseSummary(metadataForUuid.head.workflowUuid))
    metadataForUuid.foldLeft(base) {
      case (metadataSummary, metadatum) => metadataSummary |+| metadatum.toSummary
    }
  }
}

trait MetadataDatabaseAccess {
  this: ServicesStore =>

  // Only for tests
  private[impl] def closeDatabaseInterface() = databaseInterface.close()

  def addMetadataEvents(metadataEvents: Iterable[MetadataEvent])(implicit ec: ExecutionContext): Future[Unit] = {
    val metadata = metadataEvents map { metadataEvent =>
      val key = metadataEvent.key
      val workflowUuid = key.workflowId.id.toString
      val timestamp = metadataEvent.offsetDateTime.toSystemTimestamp
      val value = metadataEvent.value map { _.value }
      val valueType = metadataEvent.value map { _.valueType.typeName }
      val jobKey = key.jobKey map { jk => (jk.callFqn, jk.index, jk.attempt) }
      Metadatum(workflowUuid, key.key, jobKey.map(_._1), jobKey.flatMap(_._2), jobKey.map(_._3), value, valueType, timestamp)
    }

    databaseInterface.addMetadata(metadata)
  }

  private def metadataToMetadataEvents(workflowId: WorkflowId)(metadata: Seq[Metadatum]): Seq[MetadataEvent] = {
    metadata map { m =>
      // If callFqn is non-null then attempt will also be non-null and there is a MetadataJobKey.
      val metadataJobKey: Option[MetadataJobKey] = for {
        callFqn <- m.callFqn
        attempt <- m.attempt
      } yield MetadataJobKey(callFqn, m.index, attempt)

      val key = MetadataKey(workflowId, metadataJobKey, m.key)
      val value =  for {
        mValue <- m.value
        mType <- m.valueType
      } yield MetadataValue(mValue, MetadataType.fromString(mType))

      MetadataEvent(key, value, m.timestamp.toSystemOffsetDateTime)
    }
  }

  def queryMetadataEvents(query: MetadataQuery)(implicit ec: ExecutionContext): Future[Seq[MetadataEvent]] = {
    val uuid = query.workflowId.id.toString

    val futureMetadata: Future[Seq[Metadatum]] = query match {
      case MetadataQuery(_, None, None, None, None) => databaseInterface.queryMetadataEvents(uuid)
      case MetadataQuery(_, None, Some(key), None, None) => databaseInterface.queryMetadataEvents(uuid, key)
      case MetadataQuery(_, Some(jobKey), None, None, None) =>
        databaseInterface.queryMetadataEvents(uuid, jobKey.callFqn, jobKey.index, jobKey.attempt)
      case MetadataQuery(_, Some(jobKey), Some(key), None, None) =>
        databaseInterface.queryMetadataEvents(uuid, key, jobKey.callFqn, jobKey.index, jobKey.attempt)
      case MetadataQuery(_, None, None, Some(includeKeys), None) =>
        databaseInterface.queryMetadataEventsWithWildcardKeys(uuid, includeKeys.map(_ + "%"), requireEmptyJobKey = false)
      case MetadataQuery(_, None, None, None, Some(excludeKeys)) =>
        databaseInterface.queryMetadataEventsWithoutWildcardKeys(uuid, excludeKeys.map(_ + "%"), requireEmptyJobKey = false)
      case MetadataQuery(_, None, None, Some(includeKeys), Some(excludeKeys)) => Future.failed(
        new IllegalArgumentException(
          s"Include/Exclude keys may not be mixed: include = $includeKeys, exclude = $excludeKeys"))
      case invalidQuery => Future.failed(new IllegalArgumentException(
        s"Include/Exclude keys are only supported when querying the workflow, not when querying calls: $invalidQuery"))
    }

    futureMetadata map metadataToMetadataEvents(query.workflowId)
  }

  def queryWorkflowOutputs(id: WorkflowId)
                          (implicit ec: ExecutionContext): Future[Seq[MetadataEvent]] = {
    val uuid = id.id.toString
    databaseInterface.queryMetadataEventsWithWildcardKeys(uuid, s"${WorkflowMetadataKeys.Outputs}:%".wrapNel, requireEmptyJobKey = true).
      map(metadataToMetadataEvents(id))
  }

  def queryLogs(id: WorkflowId)
               (implicit ec: ExecutionContext): Future[Seq[MetadataEvent]] = {
    import cromwell.services.metadata.CallMetadataKeys._

    val keys = NonEmptyList(Stdout, Stderr, BackendLogsPrefix + ":%")
    databaseInterface.queryMetadataEventsWithWildcardKeys(id.id.toString, keys, requireEmptyJobKey = false) map metadataToMetadataEvents(id)
  }

  def refreshWorkflowMetadataSummaries()(implicit ec: ExecutionContext): Future[Long] = {
    databaseInterface.refreshMetadataSummaries(WorkflowMetadataKeys.StartTime, WorkflowMetadataKeys.EndTime, WorkflowMetadataKeys.Name,
      WorkflowMetadataKeys.Status, MetadataDatabaseAccess.buildUpdatedSummary)
  }

  def getWorkflowStatus(id: WorkflowId)
                       (implicit ec: ExecutionContext): Future[Option[WorkflowState]] = {
    databaseInterface.getStatus(id.toString) map { _ map WorkflowState.fromString }
  }

  def workflowExistsWithId(possibleWorkflowId: String)(implicit ec: ExecutionContext): Future[Boolean] = {
    databaseInterface.workflowExists(possibleWorkflowId)
  }

  def queryWorkflowSummaries(queryParameters: WorkflowQueryParameters)
                            (implicit ec: ExecutionContext): Future[(WorkflowQueryResponse, Option[QueryMetadata])] = {
    val workflowSummaries = databaseInterface.queryWorkflowSummaries(
      queryParameters.statuses, queryParameters.names, queryParameters.ids.map(_.toString),
      queryParameters.startDate.map(_.toSystemTimestamp), queryParameters.endDate.map(_.toSystemTimestamp),
      queryParameters.page, queryParameters.pageSize)

    val workflowSummaryCount = databaseInterface.countWorkflowSummaries(
      queryParameters.statuses, queryParameters.names, queryParameters.ids.map(_.toString),
      queryParameters.startDate.map(_.toSystemTimestamp), queryParameters.endDate.map(_.toSystemTimestamp))

    workflowSummaryCount flatMap { count =>
      workflowSummaries map { workflows =>
        (WorkflowQueryResponse(workflows.toSeq map { workflow =>
          MetadataService.WorkflowQueryResult(
            id = workflow.workflowUuid,
            name = workflow.name,
            status = workflow.status,
            start = workflow.startDate map { _.toSystemOffsetDateTime },
            end = workflow.endDate map { _.toSystemOffsetDateTime })
        }),
          //only return metadata if page is defined
          queryParameters.page map { _ => QueryMetadata(queryParameters.page, queryParameters.pageSize, Option(count)) })
      }
    }
  }
}
