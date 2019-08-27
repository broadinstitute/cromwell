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
import cromwell.services.metadata.impl.MetadataDatabaseAccess.SummaryResult
import mouse.boolean._

import scala.concurrent.duration.Duration
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
        submissionTimestamp = summary1.submissionTimestamp orElse summary2.submissionTimestamp,
        parentWorkflowExecutionUuid = summary1.parentWorkflowExecutionUuid orElse summary2.parentWorkflowExecutionUuid,
        rootWorkflowExecutionUuid = summary1.rootWorkflowExecutionUuid orElse summary2.rootWorkflowExecutionUuid,
      )
    }
  }

  def baseSummary(workflowUuid: String) =
    WorkflowMetadataSummaryEntry(workflowUuid, None, None, None, None, None, None, None)

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
        case WorkflowMetadataKeys.ParentWorkflowId =>
          base.copy(parentWorkflowExecutionUuid = metadatum.metadataValue.toRawStringOption)
        case WorkflowMetadataKeys.RootWorkflowId =>
          base.copy(rootWorkflowExecutionUuid = metadatum.metadataValue.toRawStringOption)
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

  case class SummaryResult(rowsProcessedIncreasing: Long, increasingGap: Long, rowsProcessedDecreasing: Long, decreasingGap: Long)
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

  def queryMetadataEvents(query: MetadataQuery, timeout: Duration)(implicit ec: ExecutionContext): Future[Seq[MetadataEvent]] = {

    def listKeyRequirements(keyRequirementsInput: Option[NonEmptyList[String]]): List[String] = keyRequirementsInput.map(_.toList).toList.flatten.map(_ + "%")

    val uuid = query.workflowId.id.toString

    val futureMetadata: Future[Seq[MetadataEntry]] = query match {
      case MetadataQuery(_, None, None, None, None, _) =>
        metadataDatabaseInterface.queryMetadataEntries(uuid, timeout)
      case MetadataQuery(_, None, Some(key), None, None, _) =>
        metadataDatabaseInterface.queryMetadataEntries(uuid, key, timeout)
      case MetadataQuery(_, Some(jobKey), None, None, None, _) =>
        metadataDatabaseInterface.queryMetadataEntries(uuid, jobKey.callFqn, jobKey.index, jobKey.attempt, timeout)
      case MetadataQuery(_, Some(jobKey), Some(key), None, None, _) =>
        metadataDatabaseInterface.queryMetadataEntries(uuid, key, jobKey.callFqn, jobKey.index, jobKey.attempt, timeout)
      case MetadataQuery(_, None, None, includeKeys, excludeKeys, _) =>
        val excludeKeyRequirements = listKeyRequirements(excludeKeys)
        val queryType = if (excludeKeyRequirements.contains("calls%")) WorkflowQuery else CallOrWorkflowQuery

        metadataDatabaseInterface.queryMetadataEntryWithKeyConstraints(uuid, listKeyRequirements(includeKeys), excludeKeyRequirements, queryType, timeout)
      case MetadataQuery(_, Some(MetadataQueryJobKey(callFqn, index, attempt)), None, includeKeys, excludeKeys, _) =>
        metadataDatabaseInterface.queryMetadataEntryWithKeyConstraints(uuid, listKeyRequirements(includeKeys), listKeyRequirements(excludeKeys), CallQuery(callFqn, index, attempt), timeout)
      case _ => Future.failed(new IllegalArgumentException(s"Invalid MetadataQuery: $query"))
    }

    futureMetadata map metadataToMetadataEvents(query.workflowId)
  }

  def queryWorkflowOutputs(id: WorkflowId,
                           timeout: Duration)
                          (implicit ec: ExecutionContext): Future[Seq[MetadataEvent]] = {
    val uuid = id.id.toString
    metadataDatabaseInterface.queryMetadataEntryWithKeyConstraints(
      uuid, List(s"${WorkflowMetadataKeys.Outputs}:%"), List.empty, WorkflowQuery, timeout).
      map(metadataToMetadataEvents(id))
  }

  def queryLogs(id: WorkflowId,
                timeout: Duration)
               (implicit ec: ExecutionContext): Future[Seq[MetadataEvent]] = {
    import cromwell.services.metadata.CallMetadataKeys._

    val keys = List(Stdout, Stderr, BackendLogsPrefix + ":%")
    metadataDatabaseInterface.queryMetadataEntryWithKeyConstraints(id.id.toString, keys, List.empty, CallOrWorkflowQuery, timeout) map
      metadataToMetadataEvents(id)
  }

  def refreshWorkflowMetadataSummaries(limit: Int)(implicit ec: ExecutionContext): Future[SummaryResult] = {
    for {
      (increasingProcessed, increasingGap) <- metadataDatabaseInterface.summarizeIncreasing(
        summaryNameIncreasing = WorkflowMetadataKeys.SummaryNameIncreasing,
        startMetadataKey = WorkflowMetadataKeys.StartTime,
        endMetadataKey = WorkflowMetadataKeys.EndTime,
        nameMetadataKey = WorkflowMetadataKeys.Name,
        statusMetadataKey = WorkflowMetadataKeys.Status,
        submissionMetadataKey = WorkflowMetadataKeys.SubmissionTime,
        parentWorkflowIdKey = WorkflowMetadataKeys.ParentWorkflowId,
        rootWorkflowIdKey = WorkflowMetadataKeys.RootWorkflowId,
        labelMetadataKey = WorkflowMetadataKeys.Labels,
        limit = limit,
        buildUpdatedSummary = MetadataDatabaseAccess.buildUpdatedSummary)
      (decreasingProcessed, decreasingGap) <- metadataDatabaseInterface.summarizeDecreasing(
        summaryNameDecreasing = WorkflowMetadataKeys.SummaryNameDecreasing,
        summaryNameIncreasing = WorkflowMetadataKeys.SummaryNameIncreasing,
        startMetadataKey = WorkflowMetadataKeys.StartTime,
        endMetadataKey = WorkflowMetadataKeys.EndTime,
        nameMetadataKey = WorkflowMetadataKeys.Name,
        statusMetadataKey = WorkflowMetadataKeys.Status,
        submissionMetadataKey = WorkflowMetadataKeys.SubmissionTime,
        parentWorkflowIdKey = WorkflowMetadataKeys.ParentWorkflowId,
        rootWorkflowIdKey = WorkflowMetadataKeys.RootWorkflowId,
        labelMetadataKey = WorkflowMetadataKeys.Labels,
        limit = limit,
        buildUpdatedSummary = MetadataDatabaseAccess.buildUpdatedSummary)
    } yield SummaryResult(increasingProcessed, increasingGap, decreasingProcessed, decreasingGap)
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

      val workflowLabels: Future[Map[String, String]] =
        queryParameters.additionalQueryResultFields.contains(WorkflowMetadataKeys.Labels).fold(
          metadataDatabaseInterface.getWorkflowLabels(workflow.workflowExecutionUuid), Future.successful(Map.empty))

      workflowLabels map { labels =>
        MetadataService.WorkflowQueryResult(
          id = workflow.workflowExecutionUuid,
          name = workflow.workflowName,
          status = workflow.workflowStatus,
          submission = workflow.submissionTimestamp map {_.toSystemOffsetDateTime},
          start = workflow.startTimestamp map { _.toSystemOffsetDateTime },
          end = workflow.endTimestamp map { _.toSystemOffsetDateTime },
          labels = labels.nonEmpty.option(labels),
          parentWorkflowId = workflow.parentWorkflowExecutionUuid,
          rootWorkflowId = workflow.rootWorkflowExecutionUuid
        )
      }
    }

    for {
      count <- workflowSummaryCount
      workflows <- workflowSummaries
      queryResults <- summariesToQueryResults(workflows)
    } yield (WorkflowQueryResponse(queryResults, count), queryMetadata(count))

  }
}
