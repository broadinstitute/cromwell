package cromwell.database.sql

import java.sql.Timestamp

import cats.data.NonEmptyList
import cromwell.database.sql.tables.{MetadataEntry, WorkflowMetadataSummaryEntry}

import scala.concurrent.{ExecutionContext, Future}

trait MetadataSqlDatabase {
  this: SqlDatabase =>

  /*
  The following section relates to:
.___  ___.  _______ .___________.    ___       _______       ___   .___________.    ___
|   \/   | |   ____||           |   /   \     |       \     /   \  |           |   /   \
|  \  /  | |  |__   `---|  |----`  /  ^  \    |  .--.  |   /  ^  \ `---|  |----`  /  ^  \
|  |\/|  | |   __|      |  |      /  /_\  \   |  |  |  |  /  /_\  \    |  |      /  /_\  \
|  |  |  | |  |____     |  |     /  _____  \  |  '--'  | /  _____  \   |  |     /  _____  \
|__|  |__| |_______|    |__|    /__/     \__\ |_______/ /__/     \__\  |__|    /__/     \__\
   */

  /**
    * Add metadata events to the database transactionally.
    */
  def addMetadataEntries(metadataEntries: Iterable[MetadataEntry])(implicit ec: ExecutionContext): Future[Unit]

  def metadataEntryExists(workflowExecutionUuid: String)(implicit ec: ExecutionContext): Future[Boolean]

  def queryMetadataEntries(workflowExecutionUuid: String)
                          (implicit ec: ExecutionContext): Future[Seq[MetadataEntry]]

  def queryMetadataEntries(workflowExecutionUuid: String,
                           metadataKey: String)
                          (implicit ec: ExecutionContext): Future[Seq[MetadataEntry]]

  def queryMetadataEntries(workflowExecutionUuid: String,
                           callFullyQualifiedName: String,
                           jobIndex: Option[Int],
                           jobAttempt: Int)
                          (implicit ec: ExecutionContext): Future[Seq[MetadataEntry]]

  def queryMetadataEntries(workflowUuid: String,
                           metadataKey: String,
                           callFullyQualifiedName: String,
                           jobIndex: Option[Int],
                           jobAttempt: Int)
                          (implicit ec: ExecutionContext): Future[Seq[MetadataEntry]]

  def queryMetadataEntriesLikeMetadataKeys(workflowExecutionUuid: String,
                                           metadataKeys: NonEmptyList[String],
                                           requireEmptyJobKey: Boolean)
                                          (implicit ec: ExecutionContext): Future[Seq[MetadataEntry]]

  def queryMetadataEntryNotLikeMetadataKeys(workflowExecutionUuid: String,
                                            metadataKeys: NonEmptyList[String],
                                            requireEmptyJobKey: Boolean)
                                           (implicit ec: ExecutionContext): Future[Seq[MetadataEntry]]

  /**
    * Retrieves all summarizable metadata satisfying the specified criteria.
    *
    * @param buildUpdatedSummary Takes in the optional existing summary and the metadata, returns the new summary.
    * @return A `Future` with the maximum metadataEntryId summarized by the invocation of this method.
    */
  def refreshMetadataSummaryEntries(metadataKey1: String, metadataKey2: String, metadataKey3: String, metadataKey4: String,
                                    buildUpdatedSummary:
                                    (Option[WorkflowMetadataSummaryEntry], Seq[MetadataEntry])
                                      => WorkflowMetadataSummaryEntry)
                                   (implicit ec: ExecutionContext): Future[Long]

  def getWorkflowStatus(workflowExecutionUuid: String)(implicit ec: ExecutionContext): Future[Option[String]]

  def queryWorkflowSummaries(workflowStatuses: Set[String], workflowNames: Set[String],
                             workflowExecutionUuids: Set[String], startTimestampOption: Option[Timestamp],
                             endTimestampOption: Option[Timestamp], page: Option[Int], pageSize: Option[Int])
                            (implicit ec: ExecutionContext): Future[Traversable[WorkflowMetadataSummaryEntry]]

  def countWorkflowSummaries(workflowStatuses: Set[String], workflowNames: Set[String],
                             workflowExecutionUuids: Set[String], startTimestampOption: Option[Timestamp],
                             endTimestampOption: Option[Timestamp])
                            (implicit ec: ExecutionContext): Future[Int]
}
