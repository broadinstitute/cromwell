package cromwell.database.sql

import java.sql.Timestamp

import cromwell.database.sql.tables.{Metadatum, WorkflowMetadataSummary}

import scala.concurrent.{ExecutionContext, Future}
import scalaz.NonEmptyList

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
    * Add metadata events to the database transactionally. normalized type structure is as follows:
    * (WorkflowId, MetadataKey, Option[CallFqn, CallIndex, CallAttempt], MetadataValue, MetadataValueType, Timestamp)
    */
  def addMetadata(events: Iterable[Metadatum])
                 (implicit ec: ExecutionContext): Future[Unit]

  def queryMetadataEvents(workflowUuid: String)
                         (implicit ec: ExecutionContext): Future[Seq[Metadatum]]

  def queryMetadataEvents(workflowUuid: String,
                          key: String)
                         (implicit ec: ExecutionContext): Future[Seq[Metadatum]]

  def queryMetadataEvents(workflowUuid: String,
                          callFqn: String,
                          index: Option[Int],
                          attempt: Int)
                         (implicit ec: ExecutionContext): Future[Seq[Metadatum]]

  def queryMetadataEvents(workflowUuid: String,
                          key: String,
                          callFqn: String,
                          index: Option[Int],
                          attempt: Int)
                         (implicit ec: ExecutionContext): Future[Seq[Metadatum]]

  def queryMetadataEventsWithWildcardKeys(workflowUuid: String,
                                          wildcardKeys: NonEmptyList[String],
                                          requireEmptyJobKey: Boolean)
                                         (implicit ec: ExecutionContext): Future[Seq[Metadatum]]

  def queryMetadataEventsWithoutWildcardKeys(workflowUuid: String,
                                             wildcardKeys: NonEmptyList[String],
                                             requireEmptyJobKey: Boolean)
                                            (implicit ec: ExecutionContext): Future[Seq[Metadatum]]

  /**
    * Retrieves all summarizable metadata satisfying the specified criteria.
    *
    * @param startMetadataId        The minimum ID an entry in `METADATA_JOURNAL` must have to be examined for summary.
    * @param startMetadataTimestamp An optional timestamp.  If specified, a metadatum must have a timestamp greater than
    *                               or equal to this value.
    * @param buildUpdatedSummary    Takes in the optional existing summary and the metadata, returns the new summary.
    * @return A `Future` with the metadata summarized by the invocation of this method.
    */
  def refreshMetadataSummaries(startMetadataId: Long, startMetadataTimestamp: Option[Timestamp],
                               key1: String, key2: String, key3: String, key4: String,
                               buildUpdatedSummary:
                               (Option[WorkflowMetadataSummary], Seq[Metadatum]) => WorkflowMetadataSummary)
                              (implicit ec: ExecutionContext): Future[Seq[Metadatum]]

  def getStatus(workflowUuid: String)
               (implicit ec: ExecutionContext): Future[Option[String]]

  def workflowExists(possibleWorkflowId: String)
                    (implicit ec: ExecutionContext): Future[Boolean]

  def queryWorkflowSummaries(statuses: Set[String], names: Set[String], uuids: Set[String],
                             startDate: Option[Timestamp], endDate: Option[Timestamp],
                             page: Option[Int], pageSize: Option[Int])
                            (implicit ec: ExecutionContext): Future[Traversable[WorkflowMetadataSummary]]

  def countWorkflowSummaries(statuses: Set[String], names: Set[String], uuids: Set[String],
                             startDate: Option[Timestamp], endDate: Option[Timestamp])
                            (implicit ec: ExecutionContext): Future[Int]
}
