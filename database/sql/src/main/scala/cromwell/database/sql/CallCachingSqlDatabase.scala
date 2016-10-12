package cromwell.database.sql

import cats.data.NonEmptyList
import cromwell.database.sql.joins.CallCachingJoin

import scala.concurrent.{ExecutionContext, Future}

trait CallCachingSqlDatabase {
  def addCallCaching(callCachingJoin: CallCachingJoin)(implicit ec: ExecutionContext): Future[Unit]

  def queryCallCachingEntryIds(hashKeyHashValues: NonEmptyList[(String, String)])
                              (implicit ec: ExecutionContext): Future[Seq[Int]]

  def queryCallCaching(callCachingEntryId: Int)
                      (implicit ec: ExecutionContext): Future[Option[CallCachingJoin]]

  def invalidateCall(callCachingEntryId: Int)
                    (implicit ec: ExecutionContext): Future[Unit]
}
