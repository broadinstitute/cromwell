package cromwell.engine.workflow.lifecycle.execution.callcaching

import java.nio.file.Path

import cats.data.NonEmptyList
import cromwell.backend.BackendJobExecutionActor.JobSucceededResponse
import cromwell.core.ExecutionIndex.IndexEnhancedIndex
import cromwell.core.WorkflowId
import cromwell.core.callcaching.HashResult
import cromwell.core.simpleton.WdlValueSimpleton
import cromwell.core.path.PathImplicits._
import cromwell.database.sql._
import cromwell.database.sql.joins.CallCachingJoin
import cromwell.database.sql.tables.{CallCachingDetritusEntry, CallCachingEntry, CallCachingHashEntry, CallCachingSimpletonEntry}
import cromwell.engine.workflow.lifecycle.execution.callcaching.EngineJobHashingActor.CallCacheHashes

import scala.concurrent.{ExecutionContext, Future}

final case class CallCachingEntryId(id: Int)

/**
  * Given a database-layer CallCacheStore, this accessor can access the database with engine-friendly data types.
  */
class CallCache(database: CallCachingSqlDatabase) {
  def addToCache(workflowId: WorkflowId, callCacheHashes: CallCacheHashes, response: JobSucceededResponse)(implicit ec: ExecutionContext): Future[Unit] = {
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

  private def addToCache(callCachingEntry: CallCachingEntry, hashes: Set[HashResult],
                         result: Iterable[WdlValueSimpleton], jobDetritus: Map[String, Path])
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
        case (fileName, filePath) => CallCachingDetritusEntry(fileName, filePath.toRealString)
      }
    }

    val callCachingJoin =
      CallCachingJoin(callCachingEntry, hashesToInsert.toSeq, resultToInsert.toSeq, jobDetritusToInsert.toSeq)

    database.addCallCaching(callCachingJoin)
  }

  def callCachingEntryIdsMatchingHashes(callCacheHashes: CallCacheHashes)(implicit ec: ExecutionContext): Future[Set[CallCachingEntryId]] = {
    callCachingEntryIdsMatchingHashes(NonEmptyList.fromListUnsafe(callCacheHashes.hashes.toList))
  }

  private def callCachingEntryIdsMatchingHashes(hashKeyValuePairs: NonEmptyList[HashResult])
                                       (implicit ec: ExecutionContext): Future[Set[CallCachingEntryId]] = {
    val result = database.queryCallCachingEntryIds(hashKeyValuePairs map {
      case HashResult(hashKey, hashValue) => (hashKey.key, hashValue.value)
    })

    result.map(_.toSet.map(CallCachingEntryId))
  }

  def fetchCachedResult(callCachingEntryId: CallCachingEntryId)(implicit ec: ExecutionContext): Future[Option[CallCachingJoin]] = {
    database.queryCallCaching(callCachingEntryId.id)
  }

  def invalidate(callCachingEntryId: CallCachingEntryId)(implicit ec: ExecutionContext) = {
    database.invalidateCall(callCachingEntryId.id)
  }
}
