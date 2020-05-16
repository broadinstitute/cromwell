package cromwell.backend.standard.callcaching

import com.google.common.cache.{CacheBuilder, CacheLoader}
import cromwell.core.CacheConfig
import cromwell.services.CallCaching.CallCachingEntryId

sealed abstract class BlacklistCache(val bucketCacheConfig: CacheConfig,
                                     val hitCacheConfig: CacheConfig,
                                     val name: Option[String]) {
  val bucketCache = {
    // Queries to the bucket blacklist cache return false by default (i.e. not blacklisted).
    val falseLoader = new CacheLoader[String, java.lang.Boolean]() {
      override def load(key: String): java.lang.Boolean = false
    }

    CacheBuilder.
      newBuilder().
      concurrencyLevel(bucketCacheConfig.concurrency).
      maximumSize(bucketCacheConfig.size).
      expireAfterWrite(bucketCacheConfig.ttl.length, bucketCacheConfig.ttl.unit).
      build[String, java.lang.Boolean](falseLoader)
  }

  val hitCache = {
    // Queries to the hit blacklist cache return false by default (i.e. not blacklisted).
    val falseLoader = new CacheLoader[CallCachingEntryId, java.lang.Boolean]() {
      override def load(key: CallCachingEntryId): java.lang.Boolean = false
    }

    CacheBuilder.
      newBuilder().
      concurrencyLevel(hitCacheConfig.concurrency).
      maximumSize(hitCacheConfig.size).
      expireAfterWrite(hitCacheConfig.ttl.length, hitCacheConfig.ttl.unit).
      build[CallCachingEntryId, java.lang.Boolean](falseLoader)
  }

  def isBlacklisted(hit: CallCachingEntryId): Boolean = hitCache.get(hit)

  def isBlacklisted(bucket: String): Boolean = bucketCache.get(bucket)

  def blacklistHit(hit: CallCachingEntryId): Unit = hitCache.put(hit, true)

  def blacklistBucket(bucket: String): Unit = bucketCache.put(bucket, true)
}

class RootWorkflowBlacklistCache(bucketCacheConfig: CacheConfig, hitCacheConfig: CacheConfig) extends
  BlacklistCache(bucketCacheConfig = bucketCacheConfig, hitCacheConfig = hitCacheConfig, name = None)

class GroupBlacklistCache(bucketCacheConfig: CacheConfig, hitCacheConfig: CacheConfig, val group: String) extends
  BlacklistCache(bucketCacheConfig = bucketCacheConfig, hitCacheConfig = hitCacheConfig, name = Option(group))
