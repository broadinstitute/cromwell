package cromwell.database.sql

import cats.data.NonEmptyList
import cromwell.database.sql.joins.CallCachingJoin

import scala.concurrent.{ExecutionContext, Future}

trait CallCachingSqlDatabase {
  def addCallCaching(callCachingJoin: CallCachingJoin)(implicit ec: ExecutionContext): Future[Unit]

  def queryCallCachingEntryIds(hashKeyHashValues: NonEmptyList[(String, String)])
                              (implicit ec: ExecutionContext): Future[Seq[Int]]

  def queryCallCaching(callCachingResultMetainfoId: Int)
                      (implicit ec: ExecutionContext): Future[Option[CallCachingJoin]]

  def invalidateCall(callCachingResultMetainfoId: Int)
                    (implicit ec: ExecutionContext): Future[Unit]
}
