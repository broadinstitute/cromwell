package cromwell.database.slick

import cats.data.NonEmptyList
import cromwell.database.sql._
import cromwell.database.sql.joins.CallCachingJoin
import cromwell.database.sql.tables.CallCachingEntry

import scala.concurrent.{ExecutionContext, Future}

trait CallCachingSlickDatabase extends CallCachingSqlDatabase {
  this: SlickDatabase =>

  import dataAccess.driver.api._

  override def addCallCaching(callCachingJoin: CallCachingJoin)
                             (implicit ec: ExecutionContext): Future[Unit] = {
    val action = for {
      callCachingEntryId <- dataAccess.callCachingEntryIdsAutoInc += callCachingJoin.callCachingEntry
      _ <- dataAccess.callCachingHashEntryIdsAutoInc ++= callCachingJoin.callCachingHashEntries.
        map(_.copy(callCachingEntryId = Option(callCachingEntryId)))
      _ <- dataAccess.callCachingSimpletonEntryIdsAutoInc ++= callCachingJoin.callCachingSimpletonEntries.
        map(_.copy(callCachingEntryId = Option(callCachingEntryId)))
      _ <- dataAccess.callCachingDetritusEntryIdsAutoInc ++= callCachingJoin.callCachingDetritusEntries.
        map(_.copy(callCachingEntryId = Option(callCachingEntryId)))
    } yield ()
    runTransaction(action)
  }

  override def queryCallCachingEntryIds(hashKeyHashValues: NonEmptyList[(String, String)])
                               (implicit ec: ExecutionContext): Future[Seq[Int]] = {
    val action = dataAccess.callCachingEntryIdsForHashKeyHashValues(hashKeyHashValues).result

    runTransaction(action)
  }

  override def queryCallCaching(callCachingEntryId: Int)
                               (implicit ec: ExecutionContext): Future[Option[CallCachingJoin]] = {
    val action = for {
      callCachingEntryOption <- dataAccess.
        callCachingEntriesForId(callCachingEntryId).result.headOption
      callCachingSimpletonEntries <- dataAccess.
        callCachingSimpletonEntriesForCallCachingEntryId(callCachingEntryId).result
      callCachingDetritusEntries <- dataAccess.
        callCachingDetritusEntriesForCallCachingEntryId(callCachingEntryId).result
    } yield callCachingEntryOption.map(
      CallCachingJoin(_, Seq.empty, callCachingSimpletonEntries, callCachingDetritusEntries))

    runTransaction(action)
  }

  override def invalidateCall(callCachingEntryId: Int)
                               (implicit ec: ExecutionContext): Future[Option[CallCachingEntry]] = {
    val action = for {
      _ <- dataAccess.allowResultReuseForCallCachingEntryId(callCachingEntryId).update(false)
      callCachingEntryOption <- dataAccess.callCachingEntriesForId(callCachingEntryId).result.headOption
    } yield callCachingEntryOption
    
    runTransaction(action)
  }
}
