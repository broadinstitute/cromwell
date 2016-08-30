package cromwell.database.sql

import cromwell.database.sql.tables.{CallCachingResultMetaInfoEntry, CallCachingResultSimpletonEntry}

import scala.concurrent.{ExecutionContext, Future}

trait CallCachingStore {
  def addToCache(metaInfo: CallCachingResultMetaInfoEntry, hashes: Iterable[HashKeyAndValue], result: Iterable[ResultSimpleton])(implicit ec: ExecutionContext): Future[Unit]
  def metaInfoIdsMatchingHashes(hashKeyValuePairs: Set[HashKeyAndValue])(implicit ec: ExecutionContext): Future[Set[MetaInfoId]]
  def fetchCachedResult(metaInfoId: MetaInfoId)(implicit ec: ExecutionContext): Future[Option[CachedResult]]
}

case class MetaInfoId(id: Int)
case class HashKeyAndValue(hashKey: String, hashValue: String)
case class ResultSimpleton(simpletonKey: String, simpletonValue: String, wdlType: String)
case class CachedResult(returnCode: Option[Int], resultSimpletons: Seq[CallCachingResultSimpletonEntry])
