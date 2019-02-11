package cromwell.database.sql

import java.sql.Timestamp

import cats.data.NonEmptyList
import cromwell.database.sql.joins.MetadataJobQueryValue
import cromwell.database.sql.tables.{MetadataEntry, WorkflowMetadataSummaryEntry}

import scala.concurrent.{ExecutionContext, Future}

trait MetadataSqlDatabase extends SqlDatabase {

  /*
  The following section relates to:
.___  ___.  _______ .___________.    ___       _______       ___   .___________.    ___
|   \/   | |   ____||           |   /   \     |       \     /   \  |           |   /   \
|  \  /  | |  |__   `---|  |----`  /  ^  \    |  .--.  |   /  ^  \ `---|  |----`  /  ^  \
|  |\/|  | |   __|      |  |      /  /_\  \   |  |  |  |  /  /_\  \    |  |      /  /_\  \
|  |  |  | |  |____     |  |     /  _____  \  |  '--'  | /  _____  \   |  |     /  _____  \
|__|  |__| |_______|    |__|    /__/     \__\ |_______/ /__/     \__\  |__|    /__/     \__\
   */

  def existsMetadataEntries()(implicit ec: ExecutionContext): Future[Boolean]

  /**
    * Add metadata events to the database transactionally.
    */
  def addMetadataEntries(metadataEntries: Iterable[MetadataEntry])(implicit ec: ExecutionContext): Future[Unit]

  def metadataEntryExists(workflowExecutionUuid: String)(implicit ec: ExecutionContext): Future[Boolean]

  def metadataSummaryEntryExists(workflowExecutionUuid: String)(implicit ec: ExecutionContext): Future[Boolean]

  def queryMetadataEntries(workflowExecutionUuid: String)
                          (implicit ec: ExecutionContext): Future[Seq[MetadataEntry]]

  def queryMetadataEntries(workflowExecutionUuid: String,
                           metadataKey: String)
                          (implicit ec: ExecutionContext): Future[Seq[MetadataEntry]]

  def queryMetadataEntries(workflowExecutionUuid: String,
                           callFullyQualifiedName: String,
                           jobIndex: Option[Int],
                           jobAttempt: Option[Int])
                          (implicit ec: ExecutionContext): Future[Seq[MetadataEntry]]

  def queryMetadataEntries(workflowUuid: String,
                           metadataKey: String,
                           callFullyQualifiedName: String,
                           jobIndex: Option[Int],
                           jobAttempt: Option[Int])
                          (implicit ec: ExecutionContext): Future[Seq[MetadataEntry]]

  def queryMetadataEntriesLikeMetadataKeys(workflowExecutionUuid: String,
                                           metadataKeys: NonEmptyList[String],
                                           metadataJobQueryValue: MetadataJobQueryValue)
                                          (implicit ec: ExecutionContext): Future[Seq[MetadataEntry]]

  def queryMetadataEntryNotLikeMetadataKeys(workflowExecutionUuid: String,
                                            metadataKeys: NonEmptyList[String],
                                            metadataJobQueryValue: MetadataJobQueryValue)
                                           (implicit ec: ExecutionContext): Future[Seq[MetadataEntry]]

  /**
    * Retrieves all summarizable metadata satisfying the specified criteria.
    *
    * @param buildUpdatedSummary Takes in the optional existing summary and the metadata, returns the new summary.
    * @return A `Future` with the maximum metadataEntryId summarized by the invocation of this method.
    */
  def refreshMetadataSummaryEntries(startMetadataKey: String, endMetadataKey: String, nameMetadataKey: String,
                                    statusMetadataKey: String, labelMetadataKey: String, submissionMetadataKey: String,
                                    buildUpdatedSummary:
                                    (Option[WorkflowMetadataSummaryEntry], Seq[MetadataEntry])
                                      => WorkflowMetadataSummaryEntry)
                                   (implicit ec: ExecutionContext): Future[Long]

  def getWorkflowStatus(workflowExecutionUuid: String)(implicit ec: ExecutionContext): Future[Option[String]]

  def getWorkflowLabels(workflowExecutionUuid: String)(implicit ec: ExecutionContext): Future[Map[String, String]]

  def queryWorkflowSummaries(parentIdWorkflowMetadataKey: String,
                             workflowStatuses: Set[String],
                             workflowNames: Set[String],
                             workflowExecutionUuids: Set[String],
                             labelAndKeyLabelValues: Set[(String,String)],
                             labelOrKeyLabelValues: Set[(String,String)],
                             excludeLabelAndValues: Set[(String,String)],
                             excludeLabelOrValues: Set[(String,String)],
                             submissionTimestamp: Option[Timestamp],
                             startTimestampOption: Option[Timestamp],
                             endTimestampOption: Option[Timestamp],
                             includeSubworkflows: Boolean,
                             page: Option[Int],
                             pageSize: Option[Int])
                             (implicit ec: ExecutionContext): Future[Traversable[WorkflowMetadataSummaryEntry]]

  def countWorkflowSummaries(parentIdWorkflowMetadataKey: String,
                             workflowStatuses: Set[String], workflowNames: Set[String],
                             workflowExecutionUuids: Set[String],
                             labelAndKeyLabelValues: Set[(String,String)],
                             labelOrKeyLabelValues: Set[(String,String)],
                             excludeLabelAndValues: Set[(String,String)],
                             excludeLabelOrValues: Set[(String,String)],
                             submissionTimestamp: Option[Timestamp],
                             startTimestampOption: Option[Timestamp],
                             endTimestampOption: Option[Timestamp],
                             includeSubworkflows: Boolean)
                             (implicit ec: ExecutionContext): Future[Int]
}
