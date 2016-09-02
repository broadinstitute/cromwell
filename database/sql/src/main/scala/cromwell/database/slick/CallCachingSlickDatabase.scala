package cromwell.database.slick

import cromwell.database.sql._
import cromwell.database.sql.tables.{CallCachingHashEntry, CallCachingResultMetaInfoEntry, CallCachingResultSimpletonEntry}

import scala.concurrent.{ExecutionContext, Future}

trait CallCachingSlickDatabase extends CallCachingStore {
  this: SlickDatabase =>

  import dataAccess.driver.api._

  override def addToCache(callCachingResultMetaInfo: CallCachingResultMetaInfoEntry,
                          hashesToInsert: Int => Iterable[CallCachingHashEntry],
                          resultToInsert: Int => Iterable[CallCachingResultSimpletonEntry])
                         (implicit ec: ExecutionContext): Future[Unit] = {
    val action = for {
      callCachingResultMetaInfoId <- dataAccess.callCachingResultMetaInfoAutoInc += callCachingResultMetaInfo
      _ <- dataAccess.callCachingHashAutoInc ++= hashesToInsert(callCachingResultMetaInfoId)
      _ <- dataAccess.callCachingResultSimpletonAutoInc ++= resultToInsert(callCachingResultMetaInfoId)
    } yield ()
    runTransaction(action)
  }

  override def metaInfoIdsMatchingHashes(hashKeyValuePairs: Seq[(String, String)])
                                        (implicit ec: ExecutionContext): Future[Seq[Seq[Int]]] = {

    val actions = hashKeyValuePairs map {
      case (hashKey, hashValue) => dataAccess.resultMetaInfoIdsForHashMatch(hashKey, hashValue).result
    }
    val action = DBIO.sequence(actions)

    runTransaction(action)
  }

  override def fetchCachedResult(callCachingResultMetaInfoId: Int)(implicit ec: ExecutionContext):
  Future[(Option[CallCachingResultMetaInfoEntry], Seq[CallCachingResultSimpletonEntry])] = {
    val action = for {
      metaInfo <- dataAccess.metaInfoById(callCachingResultMetaInfoId).result.headOption
      resultSimpletons <- dataAccess.resultSimpletonsForMetaInfoId(callCachingResultMetaInfoId).result
    } yield (metaInfo, resultSimpletons)

    runTransaction(action)
  }
}
