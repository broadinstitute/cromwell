package cromwell.services.metadata.impl

import cats.data.NonEmptyList
import cats.syntax.validated._
import common.validation.ErrorOr.ErrorOr
import cromwell.core._
import cromwell.database.sql.SqlConverters._
import cromwell.database.sql.joins.{CallOrWorkflowQuery, CallQuery, WorkflowQuery}
import cromwell.database.sql.tables.MetadataEntry
import cromwell.services.MetadataServicesStore
import cromwell.services.metadata._
import org.reactivestreams.Publisher
import slick.basic.DatabasePublisher

trait StreamMetadataDatabaseAccess {
  this: MetadataServicesStore =>

  private def metadataEntryToMetadataEvent(workflowId: WorkflowId)(m: MetadataEntry): MetadataEvent = {
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

  def queryMetadataEvents(query: MetadataQuery): ErrorOr[Publisher[MetadataEvent]] = {
    val uuid = query.workflowId.id.toString

    val futureMetadata: ErrorOr[DatabasePublisher[MetadataEntry]] = query match {
      case MetadataQuery(_, None, None, None, None, _) => metadataDatabaseInterface.streamQueryMetadataEntries(uuid).validNel
      case MetadataQuery(_, None, Some(key), None, None, _) => metadataDatabaseInterface.streamQueryMetadataEntries(uuid, key).validNel
      case MetadataQuery(_, Some(jobKey), None, None, None, _) =>
        metadataDatabaseInterface.streamQueryMetadataEntries(uuid, jobKey.callFqn, jobKey.index, jobKey.attempt).validNel
      case MetadataQuery(_, Some(jobKey), Some(key), None, None, _) =>
        metadataDatabaseInterface.streamQueryMetadataEntries(uuid, key, jobKey.callFqn, jobKey.index, jobKey.attempt).validNel
      case MetadataQuery(_, None, None, Some(includeKeys), None, _) =>
        metadataDatabaseInterface.
          streamQueryMetadataEntriesLikeMetadataKeys(uuid, includeKeys.map(_ + "%"), CallOrWorkflowQuery).validNel
      case MetadataQuery(_, Some(MetadataQueryJobKey(callFqn, index, attempt)), None, Some(includeKeys), None, _) =>
        metadataDatabaseInterface.
          streamQueryMetadataEntriesLikeMetadataKeys(uuid, includeKeys.map(_ + "%"), CallQuery(callFqn, index, attempt)).validNel
      case MetadataQuery(_, None, None, None, Some(excludeKeys), _) =>
        metadataDatabaseInterface.
          streamQueryMetadataEntryNotLikeMetadataKeys(uuid, excludeKeys.map(_ + "%"), CallOrWorkflowQuery).validNel
      case MetadataQuery(_, Some(MetadataQueryJobKey(callFqn, index, attempt)), None, None, Some(excludeKeys), _) =>
        metadataDatabaseInterface.
          streamQueryMetadataEntryNotLikeMetadataKeys(uuid, excludeKeys.map(_ + "%"), CallQuery(callFqn, index, attempt)).validNel
      case MetadataQuery(_, None, None, Some(includeKeys), Some(excludeKeys), _) => s"Include/Exclude keys may not be mixed: include = $includeKeys, exclude = $excludeKeys".invalidNel
      case _ => s"Invalid MetadataQuery: $query".invalidNel
    }

    futureMetadata.map(_.mapResult(metadataEntryToMetadataEvent(query.workflowId)))
  }

  def queryWorkflowOutputs(id: WorkflowId): Publisher[MetadataEvent] = {
    val uuid = id.id.toString
    metadataDatabaseInterface.streamQueryMetadataEntriesLikeMetadataKeys(
      uuid, NonEmptyList.of(s"${WorkflowMetadataKeys.Outputs}:%"), WorkflowQuery)
      .mapResult(metadataEntryToMetadataEvent(id))
  }

  def queryLogs(id: WorkflowId): Publisher[MetadataEvent] = {
    import cromwell.services.metadata.CallMetadataKeys._

    val keys = NonEmptyList.of(Stdout, Stderr, BackendLogsPrefix + ":%")
    metadataDatabaseInterface.streamQueryMetadataEntriesLikeMetadataKeys(id.id.toString, keys, CallOrWorkflowQuery).mapResult(metadataEntryToMetadataEvent(id))
  }
}
