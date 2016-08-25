package cromwell.engine.workflow.lifecycle.execution.callcaching

import cromwell.backend.BackendJobExecutionActor.SucceededResponse
import cromwell.core.ExecutionIndex.IndexEnhancedIndex
import cromwell.core.simpleton.{WdlValueBuilder, WdlValueSimpleton}
import cromwell.core.{JobOutputs, WorkflowId}
import cromwell.database.sql._
import cromwell.database.sql.tables.{CallCachingResultMetaInfoEntry, CallCachingResultSimpletonEntry}
import cromwell.engine.workflow.lifecycle.execution.callcaching.EngineJobHashingActor.CallCacheHashes
import wdl4s.TaskOutput
import wdl4s.types._
import wdl4s.values.WdlPrimitive

import scala.concurrent.{ExecutionContext, Future}
import scala.language.postfixOps

/**
  * Given a database-layer CallCacheStore, this accesser can access the database with engine-friendly data types.
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

    database.addToCache(metaInfo, hashes, result)
  }

  private def toResultSimpletons(jobOutputs: JobOutputs): Seq[ResultSimpleton] = {
    import cromwell.core.simpleton.WdlValueSimpleton._
    jobOutputs.mapValues(_.wdlValue).simplify map {
      case WdlValueSimpleton(simpletonKey, wdlPrimitive) => ResultSimpleton(simpletonKey, wdlPrimitive.valueString, wdlPrimitive.getClass.getSimpleName)
    } toSeq
  }

  def fetchMetaInfoIdsMatchingHashes(callCacheHashes: CallCacheHashes)(implicit ec: ExecutionContext): Future[Set[MetaInfoId]] = {
    database.metaInfoIdsMatchingHashes(callCacheHashes.hashes map { hash => HashKeyAndValue(hash.hashKey.key, hash.hashValue.value) })
  }

  def fetchCachedResult(metaInfoId: MetaInfoId)(implicit ec: ExecutionContext) = database.fetchCachedResult(metaInfoId)(ec)

  private def toSimpleton(entry: CallCachingResultSimpletonEntry): WdlValueSimpleton = {
    val wdlType: WdlType = entry.wdlType match {
      case "String" => WdlStringType
      case "Int" => WdlIntegerType
      case "Float" => WdlFloatType
      case "Boolean" => WdlBooleanType
      case "File" => WdlFileType
      case _ => throw new RuntimeException(s"$entry: unrecognized WDL type: ${entry.wdlType}")
    }
    WdlValueSimpleton(entry.simpletonKey, wdlType.coerceRawValue(entry.simpletonValue).get.asInstanceOf[WdlPrimitive])
  }

  def convertToJobOutputs(cachedResult: CachedResult, taskOutputs: Seq[TaskOutput])(implicit ec: ExecutionContext): JobOutputs = {
    //don't know how much error collection this *needs*
    val simpletonEntries = cachedResult.resultSimpletons map toSimpleton
    WdlValueBuilder.toJobOutputs(taskOutputs, simpletonEntries)
  }

}
