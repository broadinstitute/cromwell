package cromwell.backend.impl.htcondor.caching.provider.mongodb

import com.mongodb.DBObject
import com.mongodb.casbah.MongoCollection
import com.mongodb.casbah.commons.{MongoDBObject, TypeImports}
import com.mongodb.util.JSON
import cromwell.backend.BackendJobExecutionActor.JobSucceededResponse
import cromwell.backend.impl.htcondor.caching.CacheActor
import cromwell.backend.impl.htcondor.caching.exception.{CachedResultAlreadyExistException, CachedResultNotFoundException}
import cromwell.backend.impl.htcondor.caching.model.CachedExecutionResult
import cromwell.backend.impl.htcondor.caching.provider.mongodb.model.{KryoSerializedObject, MongoCachedExecutionResult}
import cromwell.backend.impl.htcondor.caching.provider.mongodb.serialization.KryoSerDe

class MongoCacheActor(collection: MongoCollection,
                      override val forceRewrite: Boolean = false) extends CacheActor with KryoSerDe {

  import cromwell.backend.impl.htcondor.caching.provider.mongodb.model.MongoCachedExecutionResultProtocol._
  import spray.json._

  val ErrMsg = "Got an exception when storing execution result for hash {}."
  val HashIdentifier = "hash"
  override val tag = s"[MongoCacheActor]"

  override def readExecutionResult(hash: String): CachedExecutionResult = {
    val query = MongoDBObject(HashIdentifier -> hash)
    val result = collection.findOne(query)
    result match {
      case Some(mongoDbObject) =>
        if (forceRewrite) throwCachedResultNotFoundException(hash)
        else deserializeSucceededResponse(mongoDbObject)
      case None => throwCachedResultNotFoundException(hash)
    }
  }

  override def storeExecutionResult(cachedExecutionResult: CachedExecutionResult): Unit = {
    try {
      readExecutionResult(cachedExecutionResult.hash)
      val warnMsg = s"$tag Execution result hash {${cachedExecutionResult.hash}} is already defined in database."
      log.warning(warnMsg)
      if (forceRewrite) {
        removeExecutionResult(cachedExecutionResult.hash)
        storeExecutionResultInMongoDb(cachedExecutionResult)
      }
      else throw new CachedResultAlreadyExistException(warnMsg)
    } catch {
      case e: CachedResultNotFoundException =>
        storeExecutionResultInMongoDb(cachedExecutionResult)
      case e: CachedResultAlreadyExistException =>
        throw e
      case e: Exception =>
        log.error(e, "{} Got an unhandled exception when trying to store execution result for hash {}.", tag, cachedExecutionResult.hash)
        throw e
    }
  }

  private def deserializeSucceededResponse(mongoDbObject: TypeImports.DBObject): CachedExecutionResult = {
    val cachedResult = JsonParser(mongoDbObject.toString).convertTo[MongoCachedExecutionResult]
    val succeededResponse = deserialize(cachedResult.succeededResponse.byteArray, classOf[JobSucceededResponse])
    CachedExecutionResult(cachedResult.hash, succeededResponse)
  }

  private def removeExecutionResult(hash: String): Unit = {
    val query = MongoDBObject(HashIdentifier -> hash)
    val result = collection.remove(query)
    if (result.getN == 0) throwCachedResultNotFoundException(hash)
    log.info("{} Removed execution result for hash: {}.", tag, hash)
  }

  private def storeExecutionResultInMongoDb(cachedExecutionResult: CachedExecutionResult): Unit = {
    val cachedResult = MongoCachedExecutionResult(cachedExecutionResult.hash, KryoSerializedObject(serialize(cachedExecutionResult.succeededResponse)))
    val result = collection.insert(constructDbObjectForCachedExecutionResult(cachedResult))
    if (!result.wasAcknowledged()) {
      log.error(ErrMsg, cachedExecutionResult)
      throw new IllegalStateException(ErrMsg)
    }
  }

  private def constructDbObjectForCachedExecutionResult(cachedExecutionResult: MongoCachedExecutionResult): DBObject = {
    val resultAsJsonString = cachedExecutionResult.toJson.toString()
    JSON.parse(resultAsJsonString).asInstanceOf[DBObject]
  }

  private def throwCachedResultNotFoundException(hash: String): Nothing = {
    val warnMsg = s"$tag Execution result hash {$hash} does not exist in database."
    log.warning(warnMsg)
    throw new CachedResultNotFoundException(warnMsg)
  }
}
