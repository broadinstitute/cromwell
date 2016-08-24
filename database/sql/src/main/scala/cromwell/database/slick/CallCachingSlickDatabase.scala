package cromwell.database.slick

import cromwell.database.sql._
import cromwell.database.sql.tables.{CallCachingHashEntry, CallCachingJobDetritusEntry, CallCachingResultMetaInfoEntry, CallCachingResultSimpletonEntry}

import scala.concurrent.{ExecutionContext, Future}
import scalaz.NonEmptyList

trait CallCachingSlickDatabase extends CallCachingStore {
  this: SlickDatabase =>

  import dataAccess.driver.api._

  override def addToCache(callCachingResultMetaInfo: CallCachingResultMetaInfoEntry,
                          hashesToInsert: Int => Iterable[CallCachingHashEntry],
                          resultToInsert: Int => Iterable[CallCachingResultSimpletonEntry],
                          jobDetritusToInsert: Int => Iterable[CallCachingJobDetritusEntry])
                         (implicit ec: ExecutionContext): Future[Unit] = {
    val action = for {
      callCachingResultMetaInfoId <- dataAccess.callCachingResultMetaInfoIdsAutoInc += callCachingResultMetaInfo
      _ <- dataAccess.callCachingHashEntryIdsAutoInc ++= hashesToInsert(callCachingResultMetaInfoId)
      _ <- dataAccess.callCachingResultSimpletonAutoInc ++= resultToInsert(callCachingResultMetaInfoId)
      _ <- dataAccess.callCachingJobDetritusIdAutoInc ++= jobDetritusToInsert(callCachingResultMetaInfoId)
    } yield ()
    runTransaction(action)
  }

  override def metaInfoIdsMatchingHashes(hashKeyHashValues: Seq[(String, String)])
                                        (implicit ec: ExecutionContext): Future[Seq[Seq[Int]]] = {
    val nel: NonEmptyList[(String, String)] = NonEmptyList(hashKeyHashValues.head, hashKeyHashValues.tail: _*)
    metaInfoIdsMatchingHashes(nel).map(Seq.apply(_))
  }

  def metaInfoIdsMatchingHashes(hashKeyHashValues: NonEmptyList[(String, String)])
                               (implicit ec: ExecutionContext): Future[Seq[Int]] = {
    val action = dataAccess.callCachingResultMetaInfoIdByHashKeyHashValues(hashKeyHashValues).result

    runTransaction(action)
  }

  override def fetchCachedResult(callCachingResultMetaInfoId: Int)(implicit ec: ExecutionContext):
                                  Future[(Option[CallCachingResultMetaInfoEntry],
                                  Seq[CallCachingResultSimpletonEntry],
                                  Seq[CallCachingJobDetritusEntry])] = {
    val action = for {
      metaInfo <- dataAccess.metaInfoById(callCachingResultMetaInfoId).result.headOption
      resultSimpletons <- dataAccess.resultSimpletonsForMetaInfoId(callCachingResultMetaInfoId).result
      jobDetritus <- dataAccess.jobDetritusForMetaInfoId(callCachingResultMetaInfoId).result
    } yield (metaInfo, resultSimpletons, jobDetritus)

    runTransaction(action)
  }
}
