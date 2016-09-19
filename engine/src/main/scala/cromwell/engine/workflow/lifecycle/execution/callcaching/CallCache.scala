package cromwell.engine.workflow.lifecycle.execution.callcaching

import cromwell.backend.BackendJobExecutionActor.SucceededResponse
import cromwell.core.ExecutionIndex.IndexEnhancedIndex
import cromwell.core.WorkflowId
import cromwell.core.callcaching.HashResult
import cromwell.core.simpleton.WdlValueSimpleton
import cromwell.database.sql._
import cromwell.database.sql.joins.CallCachingJoin
import cromwell.database.sql.tables.{CallCachingDetritusEntry, CallCachingEntry, CallCachingHashEntry, CallCachingSimpletonEntry}
import cromwell.engine.workflow.lifecycle.execution.callcaching.EngineJobHashingActor.CallCacheHashes

import scala.concurrent.{ExecutionContext, Future}
import scala.language.postfixOps
import scalaz.Scalaz._
import scalaz._

final case class MetaInfoId(id: Int)

case class CachedResult(callCachingEntry: CallCachingEntry, resultSimpletons: Seq[CallCachingSimpletonEntry],
                        jobDetritus: Seq[CallCachingDetritusEntry])

/**
  * Given a database-layer CallCacheStore, this accessor can access the database with engine-friendly data types.
  */
class CallCache(database: CallCachingSqlDatabase) {
  def addToCache(workflowId: WorkflowId, callCacheHashes: CallCacheHashes, response: SucceededResponse)(implicit ec: ExecutionContext): Future[Unit] = {
    val metaInfo = CallCachingEntry(
      workflowExecutionUuid = workflowId.toString,
      callFullyQualifiedName = response.jobKey.call.fullyQualifiedName,
      jobIndex = response.jobKey.index.fromIndex,
      returnCode = response.returnCode,
      allowResultReuse = true)
    val hashes = callCacheHashes.hashes
    import cromwell.core.simpleton.WdlValueSimpleton._
    val result = response.jobOutputs.mapValues(_.wdlValue).simplify
    val jobDetritus = response.jobDetritusFiles.getOrElse(Map.empty)

    addToCache(metaInfo, hashes, result, jobDetritus)
  }

  private def addToCache(metaInfo: CallCachingEntry, hashes: Set[HashResult],
                         result: Iterable[WdlValueSimpleton], jobDetritus: Map[String, String])
                        (implicit ec: ExecutionContext): Future[Unit] = {

    val hashesToInsert: Iterable[CallCachingHashEntry] = {
      hashes map { hash => CallCachingHashEntry(hash.hashKey.key, hash.hashValue.value) }
    }

    val resultToInsert: Iterable[CallCachingSimpletonEntry] = {
      result map {
        case WdlValueSimpleton(simpletonKey, wdlPrimitive) =>
          CallCachingSimpletonEntry(simpletonKey, wdlPrimitive.valueString, wdlPrimitive.wdlType.toWdlString)
      }
    }

    val jobDetritusToInsert: Iterable[CallCachingDetritusEntry] = {
      jobDetritus map {
        case (fileName, filePath) => CallCachingDetritusEntry(fileName, filePath)
      }
    }

    val callCachingJoin =
      CallCachingJoin(metaInfo, hashesToInsert.toSeq, resultToInsert.toSeq, jobDetritusToInsert.toSeq)

    database.addCallCaching(callCachingJoin)
  }

  def fetchMetaInfoIdsMatchingHashes(callCacheHashes: CallCacheHashes)(implicit ec: ExecutionContext): Future[Set[MetaInfoId]] = {
    metaInfoIdsMatchingHashes(callCacheHashes.hashes.toList.toNel.get)
  }

  private def metaInfoIdsMatchingHashes(hashKeyValuePairs: NonEmptyList[HashResult])
                                       (implicit ec: ExecutionContext): Future[Set[MetaInfoId]] = {
    val result = database.queryCallCachingEntryIds(hashKeyValuePairs map {
      case HashResult(hashKey, hashValue) => (hashKey.key, hashValue.value)
    })

    result.map(_.toSet.map(MetaInfoId))
  }

  def fetchCachedResult(metaInfoId: MetaInfoId)(implicit ec: ExecutionContext): Future[Option[CachedResult]] = {
    database.queryCallCaching(metaInfoId.id) map cachedResultOption
  }

  private def cachedResultOption(callCachingJoinOption: Option[CallCachingJoin]): Option[CachedResult] = {
    callCachingJoinOption map { callCachingJoin =>
      CachedResult(
        callCachingJoin.callCachingEntry,
        callCachingJoin.callCachingSimpletonEntries,
        callCachingJoin.callCachingDetritusEntries)
    }
  }
}
