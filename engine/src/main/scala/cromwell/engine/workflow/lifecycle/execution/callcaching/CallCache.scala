package cromwell.engine.workflow.lifecycle.execution.callcaching

import cromwell.backend.BackendJobExecutionActor.SucceededResponse
import cromwell.core.ExecutionIndex.IndexEnhancedIndex
import cromwell.core.simpleton.WdlValueSimpleton
import cromwell.core.{JobOutputs, WorkflowId}
import cromwell.database.sql._
import cromwell.database.sql.tables.{CallCachingHashEntry, CallCachingResultMetaInfoEntry, CallCachingResultSimpletonEntry}
import cromwell.engine.workflow.lifecycle.execution.callcaching.EngineJobHashingActor.CallCacheHashes

import scala.concurrent.{ExecutionContext, Future}
import scala.language.postfixOps

case class MetaInfoId(id: Int)

case class HashKeyAndValue(hashKey: String, hashValue: String)

case class ResultSimpleton(simpletonKey: String, simpletonValue: String, wdlType: String)

case class CachedResult(returnCode: Option[Int], resultSimpletons: Seq[CallCachingResultSimpletonEntry])

/**
  * Given a database-layer CallCacheStore, this accessor can access the database with engine-friendly data types.
  */
class CallCache(database: CallCachingStore) {
  def addToCache(workflowId: WorkflowId, callCacheHashes: CallCacheHashes, response: SucceededResponse)(implicit ec: ExecutionContext): Future[Unit] = {
    val metaInfo = CallCachingResultMetaInfoEntry(
      workflowUuid = workflowId.toString,
      callFqn = response.jobKey.call.fullyQualifiedName,
      scatterIndex = response.jobKey.index.fromIndex,
      returnCode = response.returnCode,
      allowResultReuse = true,
      callCachingResultMetaInfoEntryId = None)
    val hashes = callCacheHashes.hashes map { hash => HashKeyAndValue(hash.hashKey.key, hash.hashValue.value) }
    val result = toResultSimpletons(response.jobOutputs)

    addToCache(metaInfo, hashes, result)
  }

  private def addToCache(metaInfo: CallCachingResultMetaInfoEntry, hashes: Iterable[HashKeyAndValue],
                         result: Iterable[ResultSimpleton])(implicit ec: ExecutionContext): Future[Unit] = {

    def hashesToInsert(callCachingResultMetaInfoEntryId: Int): Iterable[CallCachingHashEntry] = {
      hashes map {
        case HashKeyAndValue(hashKey, hashValue) =>
          CallCachingHashEntry(hashKey, hashValue, callCachingResultMetaInfoEntryId)
      }
    }

    def resultToInsert(callCachingResultMetaInfoEntryId: Int): Iterable[CallCachingResultSimpletonEntry] = {
      result map {
        case ResultSimpleton(simpletonKey, simpletonValue, wdlType) =>
          CallCachingResultSimpletonEntry(simpletonKey, simpletonValue, wdlType, callCachingResultMetaInfoEntryId)
      }
    }

    database.addToCache(metaInfo, hashesToInsert, resultToInsert)
  }

  private def toResultSimpletons(jobOutputs: JobOutputs): Seq[ResultSimpleton] = {
    import cromwell.core.simpleton.WdlValueSimpleton._
    jobOutputs.mapValues(_.wdlValue).simplify map {
      case WdlValueSimpleton(simpletonKey, wdlPrimitive) => ResultSimpleton(simpletonKey, wdlPrimitive.valueString, wdlPrimitive.wdlType.toWdlString)
    } toSeq
  }

  def fetchMetaInfoIdsMatchingHashes(callCacheHashes: CallCacheHashes)(implicit ec: ExecutionContext): Future[Set[MetaInfoId]] = {
    metaInfoIdsMatchingHashes(callCacheHashes.hashes.toSeq map {
      hash => HashKeyAndValue(hash.hashKey.key, hash.hashValue.value)
    })
  }

  private def metaInfoIdsMatchingHashes(hashKeyValuePairs: Seq[HashKeyAndValue])
                                       (implicit ec: ExecutionContext): Future[Set[MetaInfoId]] = {
    val result = database.metaInfoIdsMatchingHashes(hashKeyValuePairs map {
      hashKeyValuePair => (hashKeyValuePair.hashKey, hashKeyValuePair.hashValue)
    })

    result.map(_.flatten.toSet.map(MetaInfoId))
  }

  def fetchCachedResult(metaInfoId: MetaInfoId)(implicit ec: ExecutionContext): Future[Option[CachedResult]] = {
    database.fetchCachedResult(metaInfoId.id) map (cachedResultOption _).tupled
  }

  private def cachedResultOption(callCachingResultMetaInfoEntryOption: Option[CallCachingResultMetaInfoEntry],
                                 callCachingResultSimpletonEntries: Seq[CallCachingResultSimpletonEntry]):
  Option[CachedResult] = {
    callCachingResultMetaInfoEntryOption map { callCachingResultMetaInfoEntry =>
      CachedResult(callCachingResultMetaInfoEntry.returnCode, callCachingResultSimpletonEntries)
    }
  }
}
