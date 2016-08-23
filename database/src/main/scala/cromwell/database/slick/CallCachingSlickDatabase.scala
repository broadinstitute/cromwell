package cromwell.database.slick

import cromwell.core.JobOutputs
import cromwell.database.sql._
import cromwell.database.sql.tables.{CallCachingHashEntry, CallCachingResultMetaInfoEntry, CallCachingResultSimpletonEntry}

import scala.concurrent.{ExecutionContext, Future}

trait CallCachingSlickDatabase extends CallCachingStore {
  this: SlickDatabase =>

  import dataAccess.driver.api._

  override def addToCache(metaInfo: CallCachingResultMetaInfoEntry, hashes: Iterable[HashKeyAndValue], result: Iterable[ResultSimpleton])(implicit ec: ExecutionContext): Future[Unit] = {
    val action = for {
      insertMetaInfoResult <- dataAccess.callCachingResultMetaInfoAutoInc += metaInfo
      metaInfoId = insertMetaInfoResult.callCachingResultMetaInfoEntryId.get
      hashesToInsert = hashes map { case HashKeyAndValue(hashKey, hashValue) => CallCachingHashEntry(hashKey, hashValue, metaInfoId, None) }
      resultToInsert = result map { case ResultSimpleton(simpletonKey, simpletonValue, wdltype) => CallCachingResultSimpletonEntry(simpletonKey, simpletonValue, wdltype, metaInfoId, None) }

      _ <- dataAccess.callCachingHashAutoInc ++= hashesToInsert
      _ <- dataAccess.callCachingResultSimpletonAutoInc ++= resultToInsert
    } yield ()
    runTransaction(action)
  }

  override def metaInfoIdsMatchingHashes(hashKeyValuePairs: Set[HashKeyAndValue])(implicit ec: ExecutionContext): Future[Set[MetaInfoId]] = {

    val actions = hashKeyValuePairs.toList map { case HashKeyAndValue(hashKey, hashValue) => dataAccess.resultMetaInfoIdsForHashMatch(hashKey, hashValue).result }
    val action = DBIO.sequence(actions)

    def setIntersection(current: Set[Int], next: Seq[Int]) = current.intersect(next.toSet)

    runTransaction(action).map { _.toList } map {
      case Nil => Set.empty
      case head :: Nil => head.toSet map MetaInfoId
      case head :: tail => tail.foldLeft(head.toSet)(setIntersection) map MetaInfoId
    }
  }

  override def fetchCachedResult(metaInfoId: MetaInfoId)(implicit ec: ExecutionContext): Future[Option[CachedResult]] = {
    val action = for {
      metaInfo <- dataAccess.metaInfoById(metaInfoId.id).result
      resultSimpletons <- dataAccess.resultSimpletonsForMetaInfoId(metaInfoId.id).result
    } yield cachedResultOption(metaInfo, resultSimpletons)

    runTransaction(action)
  }

  private def cachedResultOption(metaInfos: Seq[CallCachingResultMetaInfoEntry], simpletons: Seq[CallCachingResultSimpletonEntry]): Option[CachedResult] = {
    metaInfos.toList match {
      case Nil => None
      case head :: Nil => Option(CachedResult(head.returnCode, simpletons))
      case head :: tail => None // Impossible - unless our PK uniqueness is broken!
    }
  }


}

