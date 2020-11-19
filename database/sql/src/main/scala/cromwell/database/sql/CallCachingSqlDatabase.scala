package cromwell.database.sql

import cromwell.database.sql.joins.CallCachingJoin
import cromwell.database.sql.tables.CallCachingEntry

import scala.concurrent.{ExecutionContext, Future}

trait CallCachingSqlDatabase {
  def addCallCaching(joins: Seq[CallCachingJoin], batchSize: Int)(implicit ec: ExecutionContext): Future[Unit]

  def hasMatchingCallCachingEntriesForBaseAggregation(baseAggregationHash: String, callCachePathPrefixes: Option[List[String]])
                                                     (implicit ec: ExecutionContext): Future[Boolean]

  def findCacheHitForAggregation(baseAggregationHash: String, inputFilesAggregationHash: Option[String], callCachePathPrefixes: Option[List[String]], excludedIds: Set[Int])
                                (implicit ec: ExecutionContext): Future[Option[Int]]

  def queryResultsForCacheId(callCachingEntryId: Int)
                            (implicit ec: ExecutionContext): Future[Option[CallCachingJoin]]
  
  def callCacheJoinForCall(workflowExecutionUuid: String, callFqn: String, index: Int)
                          (implicit ec: ExecutionContext): Future[Option[CallCachingJoin]]

  def invalidateCall(callCachingEntryId: Int)
                    (implicit ec: ExecutionContext): Future[Option[CallCachingEntry]]

  def callCacheEntryIdsForWorkflowId(workflowExecutionUuid: String)(implicit ec: ExecutionContext): Future[Seq[Int]]
}
