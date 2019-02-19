package cromwell.services.metadata.impl

import cats.Semigroup
import cats.data.NonEmptyList
import cats.instances.future._
import cats.instances.list._
import cats.syntax.semigroup._
import cats.syntax.traverse._
import cromwell.core._
import cromwell.database.sql.SqlConverters._
import cromwell.database.sql.joins.{CallOrWorkflowQuery, CallQuery, WorkflowQuery}
import cromwell.database.sql.tables.{MetadataEntry, WorkflowMetadataSummaryEntry}
import cromwell.services.MetadataServicesStore
import cromwell.services.metadata.MetadataService.{QueryMetadata, WorkflowQueryResponse}
import cromwell.services.metadata._
import mouse.boolean._

import scala.concurrent.{ExecutionContext, Future}

object MetadataDatabaseAccess {

  private lazy val WorkflowMetadataSummarySemigroup = new Semigroup[WorkflowMetadataSummaryEntry] {
    override def combine(summary1: WorkflowMetadataSummaryEntry,
                         summary2: WorkflowMetadataSummaryEntry): WorkflowMetadataSummaryEntry = {
      // Resolve the status if both `this` and `that` have defined statuses.  This will evaluate to `None`
      // if one or both of the statuses is not defined.
      val resolvedStatus = for {
        thisStatus <- summary1.workflowStatus map WorkflowState.withName
        thatStatus <- summary2.workflowStatus map WorkflowState.withName
      } yield (thisStatus |+| thatStatus).toString

      WorkflowMetadataSummaryEntry(
        workflowExecutionUuid = summary1.workflowExecutionUuid,
        workflowName = summary1.workflowName orElse summary2.workflowName,
        workflowStatus = resolvedStatus orElse summary1.workflowStatus orElse summary2.workflowStatus,
        startTimestamp = summary1.startTimestamp orElse summary2.startTimestamp,
        endTimestamp = summary1.endTimestamp orElse summary2.endTimestamp,
        submissionTimestamp = summary1.submissionTimestamp orElse summary2.submissionTimestamp
      )
    }
  }

  def baseSummary(workflowUuid: String) = WorkflowMetadataSummaryEntry(workflowUuid, None, None, None, None, None)

  // If visibility is made `private`, there's a bogus warning about this being unused.
  implicit class MetadatumEnhancer(val metadatum: MetadataEntry) extends AnyVal {
    def toSummary: WorkflowMetadataSummaryEntry = {
      val base = baseSummary(metadatum.workflowExecutionUuid)
      metadatum.metadataKey match {
        case WorkflowMetadataKeys.Name => base.copy(workflowName = metadatum.metadataValue.toRawStringOption)
        case WorkflowMetadataKeys.Status => base.copy(workflowStatus = metadatum.metadataValue.toRawStringOption)
        case WorkflowMetadataKeys.StartTime =>
          base.copy(startTimestamp = metadatum.metadataValue.parseSystemTimestampOption)
        case WorkflowMetadataKeys.EndTime =>
          base.copy(endTimestamp = metadatum.metadataValue.parseSystemTimestampOption)
        case WorkflowMetadataKeys.SubmissionTime =>
          base.copy(submissionTimestamp = metadatum.metadataValue.parseSystemTimestampOption)
      }
    }
  }

  private def buildUpdatedSummary(existingSummary: Option[WorkflowMetadataSummaryEntry],
                                  metadataForUuid: Seq[MetadataEntry]): WorkflowMetadataSummaryEntry = {
    implicit val wmss = WorkflowMetadataSummarySemigroup

    val base = existingSummary.getOrElse(baseSummary(metadataForUuid.head.workflowExecutionUuid))
    metadataForUuid.foldLeft(base) {
      case (metadataSummary, metadatum) => metadataSummary |+| metadatum.toSummary
    }
  }
}

trait MetadataDatabaseAccess {
  this: MetadataServicesStore =>

  def addMetadataEvents(metadataEvents: Iterable[MetadataEvent])(implicit ec: ExecutionContext): Future[Unit] = {
    val metadata = metadataEvents map { metadataEvent =>
      val key = metadataEvent.key
      val workflowUuid = key.workflowId.id.toString
      val timestamp = metadataEvent.offsetDateTime.toSystemTimestamp
      val value = metadataEvent.value map { _.value }
      val valueType = metadataEvent.value map { _.valueType.typeName }
      val jobKey = key.jobKey map { jk => (jk.callFqn, jk.index, jk.attempt) }
      MetadataEntry(workflowUuid, jobKey.map(_._1), jobKey.flatMap(_._2), jobKey.map(_._3),
        key.key, value.toClobOption, valueType, timestamp)
    }
    metadataDatabaseInterface.addMetadataEntries(metadata)
  }

  private def metadataToMetadataEvents(workflowId: WorkflowId)(metadata: Seq[MetadataEntry]): Seq[MetadataEvent] = {
    metadata map { m =>
      // If callFullyQualifiedName is non-null then attempt will also be non-null and there is a MetadataJobKey.
      val metadataJobKey: Option[MetadataJobKey] = for {
        callFqn <- m.callFullyQualifiedName
        attempt <- m.jobAttempt
      } yield MetadataJobKey(callFqn, m.jobIndex, attempt)

      val key = MetadataKey(workflowId, metadataJobKey, m.metadataKey)
      val value = m.metadataValueType.map(mType =>
        MetadataValue(m.metadataValue.toRawString, MetadataType.fromString(mType))
      )

      MetadataEvent(key, value, m.metadataTimestamp.toSystemOffsetDateTime)
    }
  }

  def queryMetadataEvents(query: MetadataQuery)(implicit ec: ExecutionContext): Future[Seq[MetadataEvent]] = {
    val uuid = query.workflowId.id.toString

    val futureMetadata: Future[Seq[MetadataEntry]] = query match {
      case MetadataQuery(_, None, None, None, None, _) => metadataDatabaseInterface.queryMetadataEntries(uuid)
      case MetadataQuery(_, None, Some(key), None, None, _) => metadataDatabaseInterface.queryMetadataEntries(uuid, key)
      case MetadataQuery(_, Some(jobKey), None, None, None, _) =>
        metadataDatabaseInterface.queryMetadataEntries(uuid, jobKey.callFqn, jobKey.index, jobKey.attempt)
      case MetadataQuery(_, Some(jobKey), Some(key), None, None, _) =>
        metadataDatabaseInterface.queryMetadataEntries(uuid, key, jobKey.callFqn, jobKey.index, jobKey.attempt)
      case MetadataQuery(_, None, None, Some(includeKeys), None, _) =>
        metadataDatabaseInterface.
          queryMetadataEntriesLikeMetadataKeys(uuid, includeKeys.map(_ + "%"), CallOrWorkflowQuery)
      case MetadataQuery(_, Some(MetadataQueryJobKey(callFqn, index, attempt)), None, Some(includeKeys), None, _) =>
        metadataDatabaseInterface.
          queryMetadataEntriesLikeMetadataKeys(uuid, includeKeys.map(_ + "%"), CallQuery(callFqn, index, attempt))
      case MetadataQuery(_, None, None, None, Some(excludeKeys), _) =>
        metadataDatabaseInterface.
          queryMetadataEntryNotLikeMetadataKeys(uuid, excludeKeys.map(_ + "%"), CallOrWorkflowQuery)
      case MetadataQuery(_, Some(MetadataQueryJobKey(callFqn, index, attempt)), None, None, Some(excludeKeys), _) =>
        metadataDatabaseInterface.
          queryMetadataEntryNotLikeMetadataKeys(uuid, excludeKeys.map(_ + "%"), CallQuery(callFqn, index, attempt))
      case MetadataQuery(_, None, None, Some(includeKeys), Some(excludeKeys), _) => Future.failed(
        new IllegalArgumentException(
          s"Include/Exclude keys may not be mixed: include = $includeKeys, exclude = $excludeKeys"))
      case _ => Future.failed(new IllegalArgumentException(s"Invalid MetadataQuery: $query"))
    }

    futureMetadata map metadataToMetadataEvents(query.workflowId)
  }

  def queryWorkflowOutputs(id: WorkflowId)
                          (implicit ec: ExecutionContext): Future[Seq[MetadataEvent]] = {
    val uuid = id.id.toString
    metadataDatabaseInterface.queryMetadataEntriesLikeMetadataKeys(
      uuid, NonEmptyList.of(s"${WorkflowMetadataKeys.Outputs}:%"), WorkflowQuery).
      map(metadataToMetadataEvents(id))
  }

  def queryLogs(id: WorkflowId)
               (implicit ec: ExecutionContext): Future[Seq[MetadataEvent]] = {
    import cromwell.services.metadata.CallMetadataKeys._

    val keys = NonEmptyList.of(Stdout, Stderr, BackendLogsPrefix + ":%")
    metadataDatabaseInterface.queryMetadataEntriesLikeMetadataKeys(id.id.toString, keys, CallOrWorkflowQuery) map
      metadataToMetadataEvents(id)
  }

  def refreshWorkflowMetadataSummaries()(implicit ec: ExecutionContext): Future[Long] = {
    metadataDatabaseInterface.refreshMetadataSummaryEntries(WorkflowMetadataKeys.StartTime, WorkflowMetadataKeys.EndTime, WorkflowMetadataKeys.Name,
      WorkflowMetadataKeys.Status, WorkflowMetadataKeys.Labels, WorkflowMetadataKeys.SubmissionTime, MetadataDatabaseAccess.buildUpdatedSummary)
  }

  def getWorkflowStatus(id: WorkflowId)
                       (implicit ec: ExecutionContext): Future[Option[WorkflowState]] = {
    metadataDatabaseInterface.getWorkflowStatus(id.toString) map { _ map WorkflowState.withName }
  }


  def getWorkflowLabels(id: WorkflowId)(implicit ec: ExecutionContext): Future[Map[String, String]] = {
    metadataDatabaseInterface.getWorkflowLabels(id.toString)
  }

  def workflowWithIdExistsInMetadata(possibleWorkflowId: String)(implicit ec: ExecutionContext): Future[Boolean] = {
    metadataDatabaseInterface.metadataEntryExists(possibleWorkflowId)
  }

  def workflowWithIdExistsInMetadataSummaries(possibleWorkflowId: String)(implicit ec: ExecutionContext): Future[Boolean] = {
    metadataDatabaseInterface.metadataSummaryEntryExists(possibleWorkflowId)
  }

  def queryWorkflowSummaries(queryParameters: WorkflowQueryParameters)
                            (implicit ec: ExecutionContext): Future[(WorkflowQueryResponse, Option[QueryMetadata])] = {

    val labelsAndToQuery = queryParameters.labelsAnd.map(label => (label.key, label.value))
    val labelsOrToQuery = queryParameters.labelsOr.map(label => (label.key, label.value))

    val excludeLabelsAndToQuery = queryParameters.excludeLabelsAnd.map(label => (label.key, label.value))
    val excludeLabelsOrToQuery = queryParameters.excludeLabelsOr.map(label => (label.key, label.value))

    val workflowSummaries = metadataDatabaseInterface.queryWorkflowSummaries(
      WorkflowMetadataKeys.ParentWorkflowId,
      queryParameters.statuses,
      queryParameters.names,
      queryParameters.ids.map(_.toString),
      labelsAndToQuery,
      labelsOrToQuery,
      excludeLabelsAndToQuery,
      excludeLabelsOrToQuery,
      queryParameters.submissionTime.map(_.toSystemTimestamp),
      queryParameters.startDate.map(_.toSystemTimestamp),
      queryParameters.endDate.map(_.toSystemTimestamp),
      queryParameters.includeSubworkflows,
      queryParameters.page,
      queryParameters.pageSize
    )

    val workflowSummaryCount: Future[Int] = metadataDatabaseInterface.countWorkflowSummaries(
      WorkflowMetadataKeys.ParentWorkflowId,
      queryParameters.statuses,
      queryParameters.names,
      queryParameters.ids.map(_.toString),
      labelsAndToQuery,
      labelsOrToQuery,
      excludeLabelsAndToQuery,
      excludeLabelsOrToQuery,
      queryParameters.submissionTime.map(_.toSystemTimestamp),
      queryParameters.startDate.map(_.toSystemTimestamp),
      queryParameters.endDate.map(_.toSystemTimestamp),
      queryParameters.includeSubworkflows
    )

    def queryMetadata(count: Int): Option[QueryMetadata] = {
      (queryParameters.page, queryParameters.pageSize) match {
        case (None, None) => None
        case (Some(_), None) => None // Page without pagesize returns everything
        case (None, Some(_)) => Option(QueryMetadata(Option(1), queryParameters.pageSize, Option(count)))
        case _ => Option(QueryMetadata(queryParameters.page, queryParameters.pageSize, Option(count)))
      }
    }

    def summariesToQueryResults(workflows: Traversable[WorkflowMetadataSummaryEntry]): Future[List[MetadataService.WorkflowQueryResult]] = {
      workflows.toList.traverse(summaryToQueryResult)
    }

    def summaryToQueryResult(workflow: WorkflowMetadataSummaryEntry): Future[MetadataService.WorkflowQueryResult] = {

      def queryResult(labels: Map[String, String], parentId: Option[String]): MetadataService.WorkflowQueryResult = {
        MetadataService.WorkflowQueryResult(
          id = workflow.workflowExecutionUuid,
          name = workflow.workflowName,
          status = workflow.workflowStatus,
          submission = workflow.submissionTimestamp map {_.toSystemOffsetDateTime},
          start = workflow.startTimestamp map { _.toSystemOffsetDateTime },
          end = workflow.endTimestamp map { _.toSystemOffsetDateTime },
          labels = labels.nonEmpty.option(labels),
          parentWorkflowId = parentId
        )
      }

      def metadataEntriesToValue(entries: Seq[MetadataEntry]): Option[String] = {
        entries.headOption.flatMap(_.metadataValue.toRawStringOption)
      }

      def keyToMetadataValue(key: String): Future[Option[String]] = {
        val metadataEntries: Future[Seq[MetadataEntry]] = metadataDatabaseInterface.queryMetadataEntries(workflow.workflowExecutionUuid, key)
        metadataEntries map metadataEntriesToValue
      }

      def getWorkflowLabels: Future[Map[String, String]] = {
        queryParameters.additionalQueryResultFields.contains(WorkflowMetadataKeys.Labels).fold(
          metadataDatabaseInterface.getWorkflowLabels(workflow.workflowExecutionUuid), Future.successful(Map.empty))
      }

      val getParentWorkflowId: Future[Option[String]] = {
        (queryParameters.additionalQueryResultFields.contains(WorkflowMetadataKeys.ParentWorkflowId) || !queryParameters.includeSubworkflows)
          .fold(keyToMetadataValue(WorkflowMetadataKeys.ParentWorkflowId), Future.successful(None))
      }

      def formQueryResult(labels: Map[String, String], parentId: Option[String]): MetadataService.WorkflowQueryResult = queryResult(labels, parentId)


      for {
        labels <- getWorkflowLabels
        parentWorkflowId <- getParentWorkflowId
      } yield formQueryResult(labels, parentWorkflowId)
    }

    for {
      count <- workflowSummaryCount
      workflows <- workflowSummaries
      queryResults <- summariesToQueryResults(workflows)
    } yield (WorkflowQueryResponse(queryResults, count), queryMetadata(count))

  }
}
