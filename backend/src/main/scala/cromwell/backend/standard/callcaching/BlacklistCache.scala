package cromwell.backend.standard.callcaching

import com.google.common.cache.{CacheBuilder, CacheLoader}
import cromwell.core.CacheConfig
import cromwell.services.CallCaching.CallCachingEntryId

sealed trait BlacklistStatus
case object BadCacheResult extends BlacklistStatus
case object GoodCacheResult extends BlacklistStatus
case object UntestedCacheResult extends BlacklistStatus

sealed abstract class BlacklistCache(bucketCacheConfig: CacheConfig,
                                     hitCacheConfig: CacheConfig,
                                     val name: Option[String]) {
  val bucketCache = {
    // Queries to the bucket blacklist cache return UntestedCacheResult by default.
    val unknownLoader = new CacheLoader[String, BlacklistStatus]() {
      override def load(key: String): BlacklistStatus = UntestedCacheResult
    }

    CacheBuilder.
      newBuilder().
      concurrencyLevel(bucketCacheConfig.concurrency).
      maximumSize(bucketCacheConfig.size).
      expireAfterWrite(bucketCacheConfig.ttl.length, bucketCacheConfig.ttl.unit).
      build[String, BlacklistStatus](unknownLoader)
  }

  val hitCache = {
    // Queries to the hit blacklist cache return UntestedCacheResult by default (i.e. not blacklisted).
    val unknownLoader = new CacheLoader[CallCachingEntryId, BlacklistStatus]() {
      override def load(key: CallCachingEntryId): BlacklistStatus = UntestedCacheResult
    }

    CacheBuilder.
      newBuilder().
      concurrencyLevel(hitCacheConfig.concurrency).
      maximumSize(hitCacheConfig.size).
      expireAfterWrite(hitCacheConfig.ttl.length, hitCacheConfig.ttl.unit).
      build[CallCachingEntryId, BlacklistStatus](unknownLoader)
  }

  def getBlacklistStatus(hit: CallCachingEntryId): BlacklistStatus = hitCache.get(hit)

  def getBlacklistStatus(bucket: String): BlacklistStatus = bucketCache.get(bucket)

  def blacklist(hit: CallCachingEntryId): Unit = hitCache.put(hit, BadCacheResult)

  def blacklist(bucket: String): Unit = bucketCache.put(bucket, BadCacheResult)

  def whitelist(hit: CallCachingEntryId): Unit = hitCache.put(hit, GoodCacheResult)

  def whitelist(bucket: String): Unit = bucketCache.put(bucket, GoodCacheResult)
}

class RootWorkflowBlacklistCache(bucketCacheConfig: CacheConfig, hitCacheConfig: CacheConfig) extends
  BlacklistCache(bucketCacheConfig = bucketCacheConfig, hitCacheConfig = hitCacheConfig, name = None)

class GroupingBlacklistCache(bucketCacheConfig: CacheConfig, hitCacheConfig: CacheConfig, val group: String) extends
  BlacklistCache(bucketCacheConfig = bucketCacheConfig, hitCacheConfig = hitCacheConfig, name = Option(group))
