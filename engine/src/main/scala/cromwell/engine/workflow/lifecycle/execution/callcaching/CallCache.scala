package cromwell.engine.workflow.lifecycle.execution.callcaching

import cats.data.NonEmptyList
import cromwell.backend.BackendJobExecutionActor.JobSucceededResponse
import cromwell.core.ExecutionIndex.IndexEnhancedIndex
import cromwell.core.WorkflowId
import cromwell.core.callcaching.HashResult
import cromwell.core.path.Path
import cromwell.core.simpleton.WdlValueSimpleton
import cromwell.database.sql.SqlConverters._
import cromwell.database.sql._
import cromwell.database.sql.joins.CallCachingJoin
import cromwell.database.sql.tables.{CallCachingDetritusEntry, CallCachingEntry, CallCachingHashEntry, CallCachingSimpletonEntry, _}
import cromwell.engine.workflow.lifecycle.execution.callcaching.CallCache.CallCacheHashBundle
import cromwell.engine.workflow.lifecycle.execution.callcaching.CallCacheDiffQueryParameter.CallCacheDiffQueryCall
import cromwell.engine.workflow.lifecycle.execution.callcaching.CallCacheReadActor.AggregatedCallHashes
import cromwell.engine.workflow.lifecycle.execution.callcaching.EngineJobHashingActor.CallCacheHashes

import scala.concurrent.{ExecutionContext, Future}

final case class CallCachingEntryId(id: Int)
/**
  * Given a database-layer CallCacheStore, this accessor can access the database with engine-friendly data types.
  */
class CallCache(database: CallCachingSqlDatabase) {
  def addToCache(bundles: Seq[CallCacheHashBundle], batchSize: Int)(implicit ec: ExecutionContext): Future[Unit] = {
    val joins = bundles map { b =>
      val metaInfo = CallCachingEntry(
        workflowExecutionUuid = b.workflowId.toString,
        callFullyQualifiedName = b.response.jobKey.call.fullyQualifiedName,
        jobIndex = b.response.jobKey.index.fromIndex,
        jobAttempt = Option(b.response.jobKey.attempt),
        returnCode = b.response.returnCode,
        allowResultReuse = true)
      import cromwell.core.simpleton.WdlValueSimpleton._
      val result = b.response.jobOutputs.mapValues(_.wdlValue).simplify
      val jobDetritus = b.response.jobDetritusFiles.getOrElse(Map.empty)
      buildCallCachingJoin(metaInfo, b.callCacheHashes, result, jobDetritus)
    }

    database.addCallCaching(joins, batchSize)
  }

  private def buildCallCachingJoin(callCachingEntry: CallCachingEntry, callCacheHashes: CallCacheHashes,
                                   result: Iterable[WdlValueSimpleton], jobDetritus: Map[String, Path])
                                  (implicit ec: ExecutionContext): CallCachingJoin = {

    val hashesToInsert: Iterable[CallCachingHashEntry] = {
      callCacheHashes.hashes map { hash => CallCachingHashEntry(hash.hashKey.key, hash.hashValue.value) }
    }

    val aggregatedHashesToInsert: Option[CallCachingAggregationEntry] = {
      Option(CallCachingAggregationEntry(
        baseAggregation = callCacheHashes.aggregatedInitialHash,
        inputFilesAggregation = callCacheHashes.fileHashes.map(_.aggregatedHash)
      ))
    }

    val resultToInsert: Iterable[CallCachingSimpletonEntry] = {
      result map {
        case WdlValueSimpleton(simpletonKey, wdlPrimitive) =>
          CallCachingSimpletonEntry(simpletonKey, wdlPrimitive.valueString.toClobOption, wdlPrimitive.wdlType.toWdlString)
      }
    }

    val jobDetritusToInsert: Iterable[CallCachingDetritusEntry] = {
      jobDetritus map {
        case (fileName, filePath) => CallCachingDetritusEntry(fileName, filePath.pathAsString.toClobOption)
      }
    }

    CallCachingJoin(callCachingEntry, hashesToInsert.toSeq, aggregatedHashesToInsert, resultToInsert.toSeq, jobDetritusToInsert.toSeq)
  }

  def hasBaseAggregatedHashMatch(baseAggregatedHash: String)(implicit ec: ExecutionContext): Future[Boolean] = {
    database.hasMatchingCallCachingEntriesForBaseAggregation(baseAggregatedHash)
  }

  def hasKeyValuePairHashMatch(hashes: NonEmptyList[HashResult])(implicit ec: ExecutionContext): Future[Boolean] = {
    val hashKeyValuePairs = hashes map {
      case HashResult(hashKey, hashValue) => (hashKey.key, hashValue.value)
    }
    database.hasMatchingCallCachingEntriesForHashKeyValues(hashKeyValuePairs)
  }

  def callCachingHitForAggregatedHashes(aggregatedCallHashes: AggregatedCallHashes, hitNumber: Int)
                                       (implicit ec: ExecutionContext): Future[Option[CallCachingEntryId]] = {
    database.findCacheHitForAggregation(
      baseAggregationHash = aggregatedCallHashes.baseAggregatedHash,
      inputFilesAggregationHash = aggregatedCallHashes.inputFilesAggregatedHash,
      hitNumber).map(_ map CallCachingEntryId.apply)
  }

  def fetchCachedResult(callCachingEntryId: CallCachingEntryId)(implicit ec: ExecutionContext): Future[Option[CallCachingJoin]] = {
    database.queryCallCaching(callCachingEntryId.id)
  }

  def invalidate(callCachingEntryId: CallCachingEntryId)(implicit ec: ExecutionContext) = {
    database.invalidateCall(callCachingEntryId.id)
  }

  def callCacheDiff(callA: CallCacheDiffQueryCall, callB: CallCacheDiffQueryCall)(implicit ec: ExecutionContext) = {
    database.diffCallCacheHashes(
      callA.workflowId, callA.callFqn, callA.jobIndex.fromIndex,
      callB.workflowId, callB.callFqn, callB.jobIndex.fromIndex
    )
  }
}

object CallCache {
  case class CallCacheHashBundle(workflowId: WorkflowId, callCacheHashes: CallCacheHashes, response: JobSucceededResponse)
}
