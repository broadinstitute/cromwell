package cromwell.engine.workflow.lifecycle.execution.callcaching

import cromwell.backend.BackendJobExecutionActor.SucceededResponse
import cromwell.core.{JobOutputs, WorkflowId}
import cromwell.database.sql._
import cromwell.database.sql.tables.CallCachingResultMetaInfoEntry
import cromwell.engine.workflow.lifecycle.execution.callcaching.EngineJobHashingActor.CallCacheHashes
import cromwell.core.ExecutionIndex.IndexEnhancedIndex
import wdl4s.values._
import language.postfixOps

import scala.concurrent.{ExecutionContext, Future}

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
   jobOutputs flatMap { case (lqn, jobOutput) => toResultSimpletons(jobOutput.wdlValue, lqn) } toSeq
  }

  private def toResultSimpletons(wdlValue: WdlValue, simpletonKey: String): Seq[ResultSimpleton] = wdlValue match {
    case prim: WdlPrimitive => List(ResultSimpleton(simpletonKey, prim.valueString, wdlValue.getClass.getSimpleName))
    case WdlArray(_, arrayValue) => arrayValue.zipWithIndex flatMap { case (arrayItem, index) => toResultSimpletons(arrayItem, s"$simpletonKey[$index]") }
    case WdlMap(_, mapValue) => mapValue flatMap { case (key, value) => toResultSimpletons(value, s"$simpletonKey:${key.valueString}") } toSeq
    case wdlObject: WdlObjectLike => wdlObject.value flatMap { case (key, value) => toResultSimpletons(value, s"$simpletonKey:$key") } toSeq
  }

  def fetchMetaInfoIdsMatchingHashes(callCacheHashes: CallCacheHashes)(implicit ec: ExecutionContext): Future[Set[MetaInfoId]] = {
    database.metaInfoIdsMatchingHashes(callCacheHashes.hashes map { hash => HashKeyAndValue(hash.hashKey.key, hash.hashValue.value) })
  }

  def fetchCachedResult(metaInfoId: MetaInfoId)(implicit ec: ExecutionContext) = database.fetchCachedResult(metaInfoId)(ec)
}
