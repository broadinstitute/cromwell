package cromwell.database.sql

import cats.data.NonEmptyList
import cromwell.database.sql.joins.CallCachingJoin
import cromwell.database.sql.tables.CallCachingEntry

import scala.concurrent.{ExecutionContext, Future}

trait CallCachingSqlDatabase {
  def addCallCaching(joins: Seq[CallCachingJoin], batchSize: Int)(implicit ec: ExecutionContext): Future[Unit]

  def hasMatchingCallCachingEntriesForBaseAggregation(baseAggregationHash: String)
                                                     (implicit ec: ExecutionContext): Future[Boolean]

  def hasMatchingCallCachingEntriesForHashKeyValues(hashKeyHashValues: NonEmptyList[(String, String)])
                                                   (implicit ec: ExecutionContext): Future[Boolean]

  def findCacheHitForAggregation(baseAggregationHash: String, inputFilesAggregationHash: Option[String], hitNumber: Int)
                                (implicit ec: ExecutionContext): Future[Option[Int]]

  def queryResultsForCacheId(callCachingEntryId: Int)
                            (implicit ec: ExecutionContext): Future[Option[CallCachingJoin]]
  
  def cacheEntryExistsForCall(workflowExecutionUuid: String, callFqn: String, index: Int)
                             (implicit ec: ExecutionContext): Future[Boolean]

  def invalidateCall(callCachingEntryId: Int)
                    (implicit ec: ExecutionContext): Future[Option[CallCachingEntry]]
}
