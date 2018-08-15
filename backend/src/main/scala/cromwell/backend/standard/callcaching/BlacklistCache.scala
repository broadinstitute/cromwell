package cromwell.backend.standard.callcaching

import com.google.common.cache.{CacheBuilder, CacheLoader}
import cromwell.backend.standard.callcaching.BlacklistCache._
import cromwell.core.WorkflowId

import scala.concurrent.duration.FiniteDuration

class BlacklistCache(concurrency: Int, ttl: FiniteDuration) {
  val cache = {
    // Queries to the blacklist cache return false by default (i.e. not blacklisted).
    val falseLoader = new CacheLoader[BlacklistCacheKey, java.lang.Boolean]() {
      override def load(key: BlacklistCacheKey): java.lang.Boolean = false
    }

    CacheBuilder.
      newBuilder().
      concurrencyLevel(concurrency).
      expireAfterWrite(ttl.length, ttl.unit).
      build[BlacklistCacheKey, java.lang.Boolean](falseLoader)
  }

  def isBlacklisted(rootWorkflow: WorkflowId, bucket: String): Boolean = cache.get(BlacklistCacheKey(rootWorkflow, bucket))

  def blacklist(rootWorkflow: WorkflowId, bucket: String): Unit = cache.put(BlacklistCacheKey(rootWorkflow, bucket), true)
}


object BlacklistCache {
  case class BlacklistCacheKey(workflowId: WorkflowId, bucketName: String)
}
