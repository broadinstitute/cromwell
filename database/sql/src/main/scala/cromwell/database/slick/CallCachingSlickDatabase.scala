package cromwell.database.slick

import cats.data.NonEmptyList
import cats.instances.list._
import cats.instances.tuple._
import cats.syntax.foldable._
import cromwell.database.sql._
import cromwell.database.sql.joins.CallCachingJoin
import cromwell.database.sql.tables._

import scala.concurrent.{ExecutionContext, Future}

trait CallCachingSlickDatabase extends CallCachingSqlDatabase {
  this: SlickDatabase =>

  import dataAccess.driver.api._

  override def addCallCaching(joins: Seq[CallCachingJoin], batchSize: Int)
                             (implicit ec: ExecutionContext): Future[Unit] = {

    // Construct parallel lists of parent entries, hashes, simpletons, and detritus from `CallCachingJoin`s.
    val (entries, hashes, simpletons, detritus, aggregations) = joins.toList.foldMap { j =>
      (List(j.callCachingEntry), List(j.callCachingHashEntries), List(j.callCachingSimpletonEntries), List(j.callCachingDetritusEntries), List(j.callCachingAggregationEntry.toList)) }

    // Use the supplied `assigner` function to assign parent entry row IDs into the parallel `Seq` of children entities.
    def assignEntryIdsToChildren[C](ids: Seq[Int], groupingsOfChildren: Seq[Seq[C]], assigner: (Int, C) => C): Seq[C] = {
      (ids zip groupingsOfChildren) flatMap { case (id, children) => children.map(assigner(id, _)) }
    }

    // Batch insert entities into the appropriate `Table`.
    def batchInsert[E, T <: Table[E]](entries: Seq[E], tableQuery: TableQuery[T]): DBIO[_] = {
      DBIO.sequence(entries.grouped(batchSize).map { tableQuery ++= _ })
    }

    // Functions to assign call cache entry IDs into child hash entry, simpleton, and detritus rows.
    def hashAssigner(id: Int, hash: CallCachingHashEntry) = hash.copy(callCachingEntryId = Option(id))
    def simpletonAssigner(id: Int, simpleton: CallCachingSimpletonEntry) = simpleton.copy(callCachingEntryId = Option(id))
    def detritusAssigner(id: Int, detritus: CallCachingDetritusEntry) = detritus.copy(callCachingEntryId = Option(id))
    def aggregationAssigner(id: Int, aggregation: CallCachingAggregationEntry) = aggregation.copy(callCachingEntryId = Option(id))

    val action = for {
      entryIds <- dataAccess.callCachingEntryIdsAutoInc ++= entries

      hashEntries = assignEntryIdsToChildren(entryIds, hashes, hashAssigner)
      _ <- batchInsert(hashEntries, dataAccess.callCachingHashEntries)

      simpletonEntries = assignEntryIdsToChildren(entryIds, simpletons, simpletonAssigner)
      _ <- batchInsert(simpletonEntries, dataAccess.callCachingSimpletonEntries)

      detritusEntries = assignEntryIdsToChildren(entryIds, detritus, detritusAssigner)
      _ <- batchInsert(detritusEntries, dataAccess.callCachingDetritusEntries)

      aggregationEntries = assignEntryIdsToChildren(entryIds, aggregations, aggregationAssigner)
      _ <- batchInsert(aggregationEntries, dataAccess.callCachingAggregationEntries)
    } yield ()
    runTransaction(action)
  }

  override def hasMatchingCallCachingEntriesForBaseAggregation(baseAggregationHash: String)
                                   (implicit ec: ExecutionContext): Future[Boolean] = {
    val action = dataAccess.existsCallCachingEntriesForBaseAggregationHash(baseAggregationHash).result
    
    runTransaction(action)
  }

  override def findCacheHitForAggregation(baseAggregationHash: String, inputFilesAggregationHash: Option[String], hitNumber: Int)
                                         (implicit ec: ExecutionContext): Future[Option[Int]] = {
    val action = dataAccess.callCachingEntriesForAggregatedHashes(baseAggregationHash, inputFilesAggregationHash, hitNumber).result.headOption

    runTransaction(action)
  }
  
  override def hasMatchingCallCachingEntriesForHashKeyValues(hashKeyHashValues: NonEmptyList[(String, String)])
                               (implicit ec: ExecutionContext): Future[Boolean] = {
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
      CallCachingJoin(_, Seq.empty, None, callCachingSimpletonEntries, callCachingDetritusEntries))

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
