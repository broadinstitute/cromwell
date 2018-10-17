package cromwell.database.slick

import cats.data.NonEmptyList
import cromwell.database.sql.StreamMetadataSqlDatabase
import cromwell.database.sql.joins.{CallOrWorkflowQuery, CallQuery, MetadataJobQueryValue, WorkflowQuery}
import cromwell.database.sql.tables.MetadataEntry
import slick.basic.DatabasePublisher

trait StreamMetadataSlickDatabase extends StreamMetadataSqlDatabase { this: MetadataSlickDatabase =>
  import dataAccess.driver.api._

  override def streamQueryMetadataEntries(workflowExecutionUuid: String): DatabasePublisher[MetadataEntry] = {
    val action = dataAccess.metadataEntriesForWorkflowExecutionUuid(workflowExecutionUuid).result
    streamTransaction(action)
  }

  override def streamQueryMetadataEntries(workflowExecutionUuid: String,
                                    metadataKey: String): DatabasePublisher[MetadataEntry] = {
    val action =
      dataAccess.metadataEntriesForWorkflowExecutionUuidAndMetadataKey((workflowExecutionUuid, metadataKey)).result
    streamTransaction(action)
  }

  override def streamQueryMetadataEntries(workflowExecutionUuid: String,
                                    callFullyQualifiedName: String,
                                    jobIndex: Option[Int],
                                    jobAttempt: Option[Int]): DatabasePublisher[MetadataEntry] = {
    val action = dataAccess.
      metadataEntriesForJobKey((workflowExecutionUuid, callFullyQualifiedName, jobIndex, jobAttempt)).result
    streamTransaction(action)
  }

  override def streamQueryMetadataEntries(workflowUuid: String,
                                    metadataKey: String,
                                    callFullyQualifiedName: String,
                                    jobIndex: Option[Int],
                                    jobAttempt: Option[Int]): DatabasePublisher[MetadataEntry] = {
    val action = dataAccess.metadataEntriesForJobKeyAndMetadataKey((
      workflowUuid, metadataKey, callFullyQualifiedName, jobIndex, jobAttempt)).result
    streamTransaction(action)
  }

  override def streamQueryMetadataEntriesLikeMetadataKeys(workflowExecutionUuid: String,
                                                    metadataKeys: NonEmptyList[String],
                                                    metadataJobQueryValue: MetadataJobQueryValue): DatabasePublisher[MetadataEntry] = {
    val action = metadataJobQueryValue match {
      case CallQuery(callFqn, jobIndex, jobAttempt) =>
        dataAccess.metadataEntriesLikeMetadataKeysWithJob(workflowExecutionUuid, metadataKeys, callFqn, jobIndex, jobAttempt).result
      case WorkflowQuery => dataAccess.metadataEntriesLikeMetadataKeys(workflowExecutionUuid, metadataKeys, requireEmptyJobKey = true).result
      case CallOrWorkflowQuery => dataAccess.metadataEntriesLikeMetadataKeys(workflowExecutionUuid, metadataKeys, requireEmptyJobKey = false).result
    }
      
    streamTransaction(action)
  }

  override def streamQueryMetadataEntryNotLikeMetadataKeys(workflowExecutionUuid: String,
                                                     metadataKeys: NonEmptyList[String],
                                                     metadataJobQueryValue: MetadataJobQueryValue): DatabasePublisher[MetadataEntry] = {
    val action = metadataJobQueryValue match {
      case CallQuery(callFqn, jobIndex, jobAttempt) =>
        dataAccess.metadataEntriesNotLikeMetadataKeysWithJob(workflowExecutionUuid, metadataKeys, callFqn, jobIndex, jobAttempt).result
      case WorkflowQuery => dataAccess.metadataEntriesNotLikeMetadataKeys(workflowExecutionUuid, metadataKeys, requireEmptyJobKey = true).result
      case CallOrWorkflowQuery => dataAccess.metadataEntriesNotLikeMetadataKeys(workflowExecutionUuid, metadataKeys, requireEmptyJobKey = false).result
    }
    streamTransaction(action)
  }
}
