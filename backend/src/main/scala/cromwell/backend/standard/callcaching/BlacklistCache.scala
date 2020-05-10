package cromwell.backend.standard.callcaching

import com.google.common.cache.{CacheBuilder, CacheLoader}
import cromwell.core.CacheConfig
import cromwell.services.CallCaching.CallCachingEntryId

case class BlacklistCache(config: CacheConfig) {
  val bucketCache = {
    // Queries to the bucket blacklist cache return false by default (i.e. not blacklisted).
    val falseLoader = new CacheLoader[String, java.lang.Boolean]() {
      override def load(key: String): java.lang.Boolean = false
    }

    CacheBuilder.
      newBuilder().
      concurrencyLevel(config.concurrency).
      maximumSize(config.size).
      expireAfterWrite(config.ttl.length, config.ttl.unit).
      build[String, java.lang.Boolean](falseLoader)
  }

  val hitCache = {
    // Queries to the hit blacklist cache return false by default (i.e. not blacklisted).
    val falseLoader = new CacheLoader[CallCachingEntryId, java.lang.Boolean]() {
      override def load(key: CallCachingEntryId): java.lang.Boolean = false
    }

    CacheBuilder.
      newBuilder().
      concurrencyLevel(config.concurrency).
      maximumSize(config.size).
      expireAfterWrite(config.ttl.length, config.ttl.unit).
      build[CallCachingEntryId, java.lang.Boolean](falseLoader)
  }

  def isBlacklisted(hit: CallCachingEntryId, bucket: String): Boolean = hitCache.get(hit) || bucketCache.get(bucket)

  def blacklistBucket(bucket: String): Unit = bucketCache.put(bucket, true)

  def blacklistHit(hit: CallCachingEntryId): Unit = hitCache.put(hit, true)
}
